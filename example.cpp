// Config files
#include <ibamr/config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>
#include <iostream>
#include <fstream>

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/


void print_eul_data(Pointer<PatchHierarchy<NDIM> > hierarchy, LDataManager* l_data_manager);

// NOTE: These values are read in below from the input file. The values listed here are the default values.
static int post_struct_id = 0; // post structure ID number [NOTE: This assumes that the post structures are the first to
                               // be listed in the input file]
static double post_length = 6.0e-3;                   // post length [cm]
static double post_deflection_radius = 2.95e-3;       // maximum post deflection [cm]
static double post_rotational_frequency = 1.0 / .006; // post rotational frequency [1/s]

void
update_target_points(const Pointer<PatchHierarchy<NDIM> >& hierarchy,
                     const LDataManager* const l_data_manager,
                     const double current_time,
                     const double dt)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();
    const std ::pair<int, int>& post_lag_idxs =
        l_data_manager->getLagrangianStructureIndexRange(post_struct_id, finest_ln);
    Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);

    // BEG comment: I think this might be overkill --- I think we can get away with only updating the local nodes.
    // However, I also think there is no harm in updating the ghost nodes here as well.  So, let's just update
    // everything for now because we just want this to work, and it shouldn't be that much more expensive.
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    const std::vector<LNode*>& ghost_nodes = mesh->getGhostNodes();
    std::vector<LNode*> all_nodes = local_nodes;
    all_nodes.reserve(local_nodes.size() + ghost_nodes.size());
    all_nodes.insert(all_nodes.end(), ghost_nodes.begin(), ghost_nodes.end());
    for (auto* node : all_nodes)
    {
        auto* force_spec = node->getNodeDataItem<IBTargetPointForceSpec>();
        if (!force_spec) continue;
        const int lag_idx = node->getLagrangianIndex();
        Point& X_target = force_spec->getTargetPointPosition();
        if (post_lag_idxs.first <= lag_idx && lag_idx < post_lag_idxs.second)
        {
            // We are using the z component of the targeted position to figure out where we are on the post. This
            // fundamentally assumes that the posts are anchored at z=0.  If the posts are anchored at z != 0, then that
            // needs to be taken into account here.
            //
            // This works because we DO NOT change the value of X[2]!
            //
            // Also note that if the post is not initialized in the "leaning" configuration, then these forces will act
            // to "stretch" out the posts.  We could try to account for that here, but it is simplest to initialize the
            // post configuration (in the .vertex file) in the leaning configuration consistent with this forcing.
            double h = X_target[2];
            if (h > sqrt(std::numeric_limits<double>::epsilon()))
            {
                // The "slanted post height" is the maximum z value in the slanted configuration. It is determined from
                // the post length and the prescribed deflection radius.
                //
                // Notice that we do not pre-compute this value above because doing so would make it harder to read in
                // user defined post lengths & deflection radii.
                double slanted_post_height =
                    sqrt(post_length * post_length - post_deflection_radius * post_deflection_radius);
                double deflection = post_deflection_radius * h / slanted_post_height;
                double previous_time = current_time - dt;
                X_target[0] = X_target[0] - deflection * cos(2.0 * M_PI * previous_time * post_rotational_frequency) + deflection * cos(2.0 * M_PI * current_time * post_rotational_frequency);
                X_target[1] = X_target[1] - deflection * sin(2.0 * M_PI * previous_time * post_rotational_frequency) + deflection * sin(2.0 * M_PI * current_time * post_rotational_frequency);

            }
        }
    }
}



int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool is_from_restart = app_initializer->isFromRestart();
        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Determine whether to include an dye concentration field.
        const bool simulate_dye_concentration_field =
            input_db->getBoolWithDefault("SIMULATE_DYE_CONCENTRATION_FIELD", false);

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        if (simulate_dye_concentration_field)
        {
            const string adv_diff_solver_type = input_db->getStringWithDefault("ADV_DIFF_SOLVER_TYPE", "SEMI_IMPLICIT");
            if (adv_diff_solver_type == "GODUNOV")
            {
                Pointer<AdvectorExplicitPredictorPatchOps> predictor = new AdvectorExplicitPredictorPatchOps(
                    "AdvectorExplicitPredictorPatchOps",
                    app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
                adv_diff_integrator = new AdvDiffPredictorCorrectorHierarchyIntegrator(
                    "AdvDiffPredictorCorrectorHierarchyIntegrator",
                    app_initializer->getComponentDatabase("AdvDiffPredictorCorrectorHierarchyIntegrator"),
                    predictor);
            }
            else if (adv_diff_solver_type == "SEMI_IMPLICIT")
            {
                adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                    "AdvDiffSemiImplicitHierarchyIntegrator",
                    app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
            }
            else
            {
                TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n"
                                                       << "Valid options are: GODUNOV, SEMI_IMPLICIT");
            }
            navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        }
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = nullptr;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
        }

        if (simulate_dye_concentration_field)
        {
            // Setup the advected and diffused quantity.
            Pointer<CellVariable<NDIM, double> > Q_var = new CellVariable<NDIM, double>("Q");
            adv_diff_integrator->registerTransportedQuantity(Q_var);
            adv_diff_integrator->setDiffusionCoefficient(Q_var, input_db->getDouble("D"));
            if (input_db->keyExists("ConcentrationInitialConditions"))
            {
                Pointer<CartGridFunction> Q_init = new muParserCartGridFunction(
                    "Q_init", app_initializer->getComponentDatabase("ConcentrationInitialConditions"), grid_geometry);
                adv_diff_integrator->setInitialConditions(Q_var, Q_init);
            }
            Pointer<FaceVariable<NDIM, double> > u_adv_var = navier_stokes_integrator->getAdvectionVelocityVariable();
            adv_diff_integrator->setAdvectionVelocity(Q_var, u_adv_var);
            std::unique_ptr<RobinBcCoefStrategy<NDIM> > Q_bc_coefs = nullptr;
            if (periodic_shift.min() == 0)
            {
                Q_bc_coefs.reset(new muParserRobinBcCoefs(
                    "Q_bc_coefs", app_initializer->getComponentDatabase("ConcentrationBcCoefs"), grid_geometry));
                adv_diff_integrator->setPhysicalBcCoef(Q_var, Q_bc_coefs.get());
            }
        }

        // Get information on body forcing.
        const double rot_frequency = input_db->getDouble("ROT_FREQUENCY");
        const double post_force = input_db->getDouble("POST_FORCE");

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write restart data before starting main time integration loop.
        if (dump_restart_data && !is_from_restart)
        {
            pout << "\nWriting restart files...\n\n";
            RestartManager::getManager()->writeRestartFile(restart_dump_dirname, 0);
        }

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }

        // Create Text File which Stores the locations of each IB point

        ofstream fout("Coordinates.txt",ios::out);




        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            #if 0
                        IBTK::Vector F;
                        F(0) = post_force * cos(loop_time * rot_frequency * 2.0 * M_PI);
                        F(1) = post_force * sin(loop_time * rot_frequency * 2.0 * M_PI);
                        F(2) = 0.0;
                        pout << F << "\n";
                        ib_force_fcn->setUniformBodyForce(F, /*structure_id*/ 0, patch_hierarchy->getFinestLevelNumber());
            #else
                        update_target_points(patch_hierarchy, ib_method_ops->getLDataManager(), loop_time, dt);
            #endif

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            print_eul_data(patch_hierarchy, ib_method_ops->getLDataManager());
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            ofstream fout("Coordinates.txt",ios::app);

            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<LData> X_data = ib_method_ops->getLDataManager()->getLData("X",finest_ln);
            Vec X_vec = X_data->getVec();
            double* x_vals;
            int ierr = VecGetArray(X_vec, &x_vals);
            IBTK_CHKERRQ(ierr);
            Pointer<LMesh> l_mesh = ib_method_ops->getLDataManager()->getLMesh(finest_ln);
            const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
            for (const auto& node : local_nodes)
                {
                    const int lag_idx = node->getLagrangianIndex();
                    const int petsc_idx = node ->getLocalPETScIndex();
                    Eigen::Map<Vector3d> X(&x_vals[petsc_idx*NDIM]);
                    fout.unsetf(ios_base::showpos);
                    fout.setf(ios_base::scientific);
                    fout.precision(5);
                    fout << loop_time << "," << X(0) <<","<< X(1) << "," << X(2) <<"\n";
                }


            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
            }
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main

void print_eul_data(Pointer<PatchHierarchy<NDIM>> hierarchy, LDataManager* l_data_manager)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();
    Pointer<LData> X_data = l_data_manager->getLData("X",finest_ln);
    Vec X_vec = X_data->getVec();
    double* x_vals;
    int ierr = VecGetArray(X_vec, &x_vals);
    IBTK_CHKERRQ(ierr);
    Pointer<LMesh> l_mesh = l_data_manager->getLMesh(finest_ln);
    const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
    for (const auto& node : local_nodes)
    {
        const int lag_idx = node->getLagrangianIndex();
        const int petsc_idx = node ->getLocalPETScIndex();
        Eigen::Map<Vector3d> X(&x_vals[petsc_idx*NDIM]);
        pout << "Eulerian Location of node" << lag_idx << "\n";
        pout << X <<"\n";
    }
    return;
}
