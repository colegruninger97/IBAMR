// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/MarkerPatchHierarchy.h>

#include <CartesianGridGeometry.h>
#include <mpi.h>

#include <numeric>

#include <ibtk/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace
{
std::vector<std::vector<hier::Box<NDIM> > >
compute_nonoverlapping_patch_boxes(const Pointer<BasePatchLevel<NDIM> >& c_level,
                                   const Pointer<BasePatchLevel<NDIM> >& f_level)
{
    const Pointer<PatchLevel<NDIM> > coarse_level = c_level;
    const Pointer<PatchLevel<NDIM> > fine_level = f_level;
    TBOX_ASSERT(coarse_level);
    TBOX_ASSERT(fine_level);
    TBOX_ASSERT(coarse_level->getLevelNumber() + 1 == fine_level->getLevelNumber());

    const IntVector<NDIM> ratio = fine_level->getRatioToCoarserLevel();

    // Get all (including those not on this processor) fine-level boxes:
    BoxList<NDIM> finer_box_list;
    long combined_size = 0;
    for (int i = 0; i < fine_level->getNumberOfPatches(); ++i)
    {
        Box<NDIM> patch_box = fine_level->getBoxForPatch(i);
        patch_box.coarsen(ratio);
        combined_size += patch_box.size();
        finer_box_list.addItem(patch_box);
    }
    finer_box_list.simplifyBoxes();

    // Remove said boxes from each coarse-level patch:
    const auto rank = IBTK_MPI::getRank();
    std::vector<std::vector<Box<NDIM> > > result;
    long coarse_size = 0;
    for (int i = 0; i < coarse_level->getNumberOfPatches(); ++i)
    {
        BoxList<NDIM> coarse_box_list;
        coarse_box_list.addItem(coarse_level->getBoxForPatch(i));
        coarse_size += coarse_box_list.getFirstItem().size();
        coarse_box_list.removeIntersections(finer_box_list);

        const bool patch_is_local = rank == coarse_level->getMappingForPatch(i);
        if (patch_is_local) result.emplace_back();
        typename tbox::List<Box<NDIM> >::Iterator it(coarse_box_list);
        while (it)
        {
            if (patch_is_local) result.back().push_back(*it);
            combined_size += (*it).size();
            it++;
        }
    }

    TBOX_ASSERT(coarse_size == combined_size);

    return result;
}

void
do_interpolation(const int data_idx,
                 const std::vector<double>& positions,
                 const Pointer<Patch<NDIM> > patch,
                 const std::string& kernel,
                 std::vector<double>& velocities)
{
    Pointer<PatchData<NDIM> > data = patch->getPatchData(data_idx);
    Pointer<CellData<NDIM, double> > cc_data = data;
    Pointer<SideData<NDIM, double> > sc_data = data;
    const bool is_cc_data = cc_data;
    const bool is_sc_data = sc_data;
    // Only interpolate things within 1 cell of the patch box - we aren't
    // guaranteed to have more ghost data than that
    Box<NDIM> interp_box = data->getBox();
    interp_box.grow(1);
#ifndef NDEBUG
    std::fill(velocities.begin(), velocities.end(), std::numeric_limits<double>::signaling_NaN());
#endif

    if (is_cc_data)
    {
        LEInteractor::interpolate(velocities, NDIM, positions, NDIM, cc_data, patch, interp_box, kernel);
    }
    else if (is_sc_data)
    {
        LEInteractor::interpolate(velocities, NDIM, positions, NDIM, sc_data, patch, interp_box, kernel);
    }
    else
    {
        TBOX_ERROR("not implemented");
    }
#ifndef NDEBUG
    for (const double& v : velocities)
    {
        if (std::isnan(v))
        {
            TBOX_ERROR(
                "One of the marker points inside a patch did not have its velocity set by interpolation. "
                "This most likely means that the marker point is outside of its assigned patch, which "
                "should not happen at this point.");
        }
    }
#endif
}
} // namespace
MarkerPatch::MarkerPatch(const Box<NDIM>& patch_box,
                         const std::vector<Box<NDIM> >& nonoverlapping_patch_boxes,
                         const Pointer<CartesianGridGeometry<NDIM> >& grid_geom,
                         const IntVector<NDIM>& ratio)
    : d_patch_box(patch_box), d_nonoverlapping_patch_boxes(nonoverlapping_patch_boxes)
{
    for (const Box<NDIM>& box : nonoverlapping_patch_boxes)
    {
        TBOX_ASSERT(patch_box.contains(box));
    }
    d_domain_box = Box<NDIM>::refine(grid_geom->getPhysicalDomain()[0], ratio);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_x_lo[d] = grid_geom->getXLower()[d];
        d_x_up[d] = grid_geom->getXUpper()[d];
        d_dx[d] = grid_geom->getDx()[d] / ratio(d);
    }
}

void
MarkerPatch::insert(const int& index, const IBTK::Point& position, const IBTK::Vector& velocity)
{
    TBOX_ASSERT(index >= 0);
    const auto p = std::lower_bound(d_indices.begin(), d_indices.end(), index) - d_indices.begin();
    d_indices.insert(d_indices.begin() + p, index);
    d_positions.insert(d_positions.begin() + NDIM * p, position.data(), position.data() + position.size());
    d_velocities.insert(d_velocities.begin() + NDIM * p, velocity.data(), velocity.data() + velocity.size());
}

bool
MarkerPatch::contains(const IBTK::Point& position) const
{
    const auto index = IndexUtilities::getCellIndex(
        position.data(), d_x_lo.data(), d_x_up.data(), d_dx.data(), d_domain_box.lower(), d_domain_box.upper());
    if (!d_patch_box.contains(index)) return false;

    for (const auto& box : d_nonoverlapping_patch_boxes)
        if (box.contains(index)) return true;
    return false;
}

std::tuple<std::vector<int>, EigenAlignedVector<IBTK::Point>, EigenAlignedVector<IBTK::Vector> >
MarkerPatch::prune()
{
    std::vector<int> indices;
    EigenAlignedVector<IBTK::Point> positions;
    EigenAlignedVector<IBTK::Vector> velocities;
    // Prune from the back to the front to avoid moving stored markers
    for (int index = size() - 1; index >= 0; index--)
    {
        IBTK::Point point;
        IBTK::Point velocity;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            point[d] = d_positions[index * NDIM + d];
            velocity[d] = d_velocities[index * NDIM + d];
        }
        if (!contains(point))
        {
            // emplace so that things remain sorted
            indices.insert(indices.begin(), d_indices[index]);
            d_indices.erase(d_indices.begin() + index);
            auto p0 = d_positions.begin() + NDIM * index;
            auto p1 = d_positions.begin() + NDIM * (index + 1);
            positions.emplace(positions.begin(), &*p0);
            d_positions.erase(p0, p1);
            auto v0 = d_velocities.begin() + NDIM * index;
            auto v1 = d_velocities.begin() + NDIM * (index + 1);
            velocities.emplace(velocities.begin(), &*v0);
            d_velocities.erase(v0, v1);
        }
    }

    return std::make_tuple(std::move(indices), std::move(positions), std::move(velocities));
}

std::tuple<int, IBTK::Point, IBTK::Vector>
MarkerPatch::operator[](const unsigned int local_index) const
{
#ifndef NDEBUG
    TBOX_ASSERT(local_index < size());
#endif
    IBTK::Point position(&d_positions[local_index * NDIM]);
    IBTK::Point velocity(&d_velocities[local_index * NDIM]);
    return std::make_tuple(d_indices[local_index], position, velocity);
}

std::size_t
MarkerPatch::size() const
{
#ifndef NDEBUG
    TBOX_ASSERT(d_positions.size() == d_indices.size() * NDIM);
    TBOX_ASSERT(d_velocities.size() == d_indices.size() * NDIM);
#endif
    return d_indices.size();
}

MarkerPatchHierarchy::MarkerPatchHierarchy(const std::string& name,
                                           Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                           const EigenAlignedVector<IBTK::Point>& positions,
                                           const EigenAlignedVector<IBTK::Point>& velocities)
    : d_object_name(name),
      d_num_markers(positions.size()),
      d_hierarchy(patch_hierarchy),
      d_markers_outside_domain(
          Box<NDIM>(std::numeric_limits<int>::lowest(), std::numeric_limits<int>::max()),
          std::vector<Box<NDIM> >{ Box<NDIM>(std::numeric_limits<int>::lowest(), std::numeric_limits<int>::max()) },
          d_hierarchy->getGridGeometry(),
          IntVector<NDIM>(1))
{
    reinit(positions, velocities);
}

void
MarkerPatchHierarchy::reinit(const EigenAlignedVector<IBTK::Point>& positions,
                             const EigenAlignedVector<IBTK::Point>& velocities)
{
    TBOX_ASSERT(positions.size() == velocities.size());
    const auto rank = IBTK_MPI::getRank();
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    unsigned int num_emplaced_markers = 0;
    std::vector<bool> marker_emplaced(positions.size());

    auto insert_markers = [&](MarkerPatch& marker_patch) {
        for (unsigned int k = 0; k < positions.size(); ++k)
        {
            if (!marker_emplaced[k] && marker_patch.contains(positions[k]))
            {
                marker_emplaced[k] = true;
                marker_patch.insert(k, positions[k], velocities[k]);
                ++num_emplaced_markers;
            }
        }
    };

    d_marker_patches.clear();
    d_marker_patches.resize(d_hierarchy->getFinestLevelNumber() + 1);

    // Assign particles to levels in a top-down way:
    //
    // 1. Compute a set of boxes for the patch which do not intersect any
    //    finer patch level.
    // 2. Emplace particles in the vector corresponding to the local index of a
    //    patch in the present level.
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        Pointer<PatchLevel<NDIM> > finer_level =
            ln == d_hierarchy->getFinestLevelNumber() ? nullptr : d_hierarchy->getPatchLevel(ln + 1);
        const IntVector<NDIM>& ratio = current_level->getRatio();

        // If there is no finer level then each Patch has exactly one
        // nonoverlapping box
        if (!finer_level)
        {
            for (int i = 0; i < current_level->getNumberOfPatches(); ++i)
            {
                if (rank == current_level->getMappingForPatch(i))
                {
                    const Box<NDIM> box = current_level->getPatch(i)->getBox();
                    d_marker_patches[ln].emplace_back(box, std::vector<Box<NDIM> >{ box }, grid_geom, ratio);
                    insert_markers(d_marker_patches[ln].back());
                }
            }
        }
        // otherwise we need to subtract off the boxes on the finer level first.
        else
        {
            const std::vector<std::vector<hier::Box<NDIM> > > nonoverlapping_patch_boxes =
                compute_nonoverlapping_patch_boxes(current_level, finer_level);
            unsigned int local_num = 0;
            for (int i = 0; i < current_level->getNumberOfPatches(); ++i)
            {
                if (rank == current_level->getMappingForPatch(i))
                {
                    d_marker_patches[ln].emplace_back(
                        current_level->getPatch(i)->getBox(), nonoverlapping_patch_boxes[local_num], grid_geom, ratio);
                    insert_markers(d_marker_patches[ln].back());
                    ++local_num;
                }
            }
        }
    }

    // Handle any markers which may be outside of the domain:
    if (rank == 0)
    {
        d_markers_outside_domain.d_indices.resize(0);
        d_markers_outside_domain.d_positions.resize(0);
        d_markers_outside_domain.d_velocities.resize(0);

        const Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
        const double* const domain_x_lower = grid_geom->getXLower();
        const double* const domain_x_upper = grid_geom->getXUpper();
        for (unsigned int k = 0; k < positions.size(); ++k)
        {
            const auto& position = positions[k];
            bool point_outside_domain = false;
            for (unsigned int d = 0; d < NDIM; ++d)
                point_outside_domain =
                    point_outside_domain || ((position[d] < domain_x_lower[d]) || (domain_x_upper[d] <= position[d]));

            if (point_outside_domain)
            {
                marker_emplaced[k] = true;
                IBTK::Vector v;
                v.fill(0);
                d_markers_outside_domain.insert(k, position, v);
                ++num_emplaced_markers;
            }
        }
    }

    num_emplaced_markers = IBTK_MPI::sumReduction(num_emplaced_markers);
    TBOX_ASSERT(num_emplaced_markers == positions.size());
}

const MarkerPatch&
MarkerPatchHierarchy::getMarkerPatch(const int ln, const int local_patch_num) const
{
    TBOX_ASSERT(ln < static_cast<int>(d_marker_patches.size()));
    TBOX_ASSERT(local_patch_num < static_cast<int>(d_marker_patches[ln].size()));
    return d_marker_patches[ln][local_patch_num];
}

std::size_t
MarkerPatchHierarchy::getNumberOfMarkers() const
{
    return d_num_markers;
}

void
MarkerPatchHierarchy::setVelocities(const int u_idx, const std::string& kernel)
{
    const auto rank = IBTK_MPI::getRank();
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        unsigned int local_patch_num = 0;
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                Pointer<Patch<NDIM> > patch = current_level->getPatch(p);
                MarkerPatch& marker_patch = d_marker_patches[ln][local_patch_num];
                do_interpolation(u_idx, marker_patch.d_positions, patch, kernel, marker_patch.d_velocities);
                ++local_patch_num;
            }
        }
    }
}

void
MarkerPatchHierarchy::midpointStep(const double dt,
                                   const int u_half_idx,
                                   const int u_new_idx,
                                   const std::string& kernel)
{
    const auto rank = IBTK_MPI::getRank();
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        unsigned int local_patch_num = 0;
        Pointer<PatchLevel<NDIM> > current_level = d_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                MarkerPatch& marker_patch = d_marker_patches[ln][local_patch_num];
                Pointer<Patch<NDIM> > patch = current_level->getPatch(p);

                // 1. Do a forward Euler step:
                std::vector<double> half_positions(marker_patch.d_positions);
                for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
                {
                    half_positions[i] = marker_patch.d_positions[i] + 0.5 * dt * marker_patch.d_velocities[i];
                }

                // 2. Interpolate midpoint velocity:
                std::vector<double> half_velocities(marker_patch.d_velocities);
                do_interpolation(u_half_idx, half_positions, patch, kernel, half_velocities);

                // 3. Do a midpoint step:
                for (unsigned int i = 0; i < marker_patch.d_positions.size(); ++i)
                {
                    marker_patch.d_positions[i] += dt * half_velocities[i];
                }

                // 4. Interpolate the velocity at the new time:
                do_interpolation(u_new_idx, marker_patch.d_positions, patch, kernel, marker_patch.d_velocities);

                ++local_patch_num;
            }
        }
    }

    pruneAndRedistribute();
}

void
MarkerPatchHierarchy::pruneAndRedistribute()
{
    // 1. Collect all markers which have left their respective patches:
    std::vector<int> moved_indices;
    std::vector<double> moved_positions;
    std::vector<double> moved_velocities;
    for (auto& level_marker_patches : d_marker_patches)
    {
        for (auto& marker_patch : level_marker_patches)
        {
            const auto moved_markers = marker_patch.prune();
            const auto n_moved_markers = std::get<0>(moved_markers).size();
            moved_indices.reserve(moved_indices.size() + n_moved_markers);
            moved_positions.reserve(moved_indices.size() + n_moved_markers);
            moved_velocities.reserve(moved_indices.size() + n_moved_markers);
            for (unsigned int k = 0; k < std::get<0>(moved_markers).size(); ++k)
            {
                moved_indices.push_back(std::get<0>(moved_markers)[k]);
                moved_positions.insert(moved_positions.end(),
                                       std::get<1>(moved_markers)[k].data(),
                                       std::get<1>(moved_markers)[k].data() + NDIM);
                moved_velocities.insert(moved_velocities.end(),
                                        std::get<2>(moved_markers)[k].data(),
                                        std::get<2>(moved_markers)[k].data() + NDIM);
            }
        }
    }

    // 2. Communicate data so all processors have all errant marker points:
    const auto rank = IBTK_MPI::getRank();
    const auto n_procs = IBTK_MPI::getNodes();
    const int num_indices = moved_indices.size();
    std::vector<int> counts(n_procs);
    int ierr = MPI_Allgather(&num_indices, 1, MPI_INT, counts.data(), 1, MPI_INT, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);
    // We want the first entry to be zero and the last to be the sum
    std::vector<int> offsets(n_procs + 1);
    std::partial_sum(counts.begin(), counts.end(), offsets.begin() + 1);
    std::vector<int> vector_counts;
    std::vector<int> vector_offsets;
    for (int r = 0; r < n_procs; ++r)
    {
        vector_counts.push_back(counts[r] * NDIM);
        vector_offsets.push_back(offsets[r] * NDIM);
    }
    const auto num_total_indices = offsets.back();

    std::vector<int> new_indices(num_total_indices);
    std::vector<double> new_positions(num_total_indices * NDIM);
    std::vector<double> new_velocities(num_total_indices * NDIM);
    ierr = MPI_Allgatherv(moved_indices.data(),
                          num_indices,
                          MPI_INT,
                          new_indices.data(),
                          counts.data(),
                          offsets.data(),
                          MPI_INT,
                          IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);
    ierr = MPI_Allgatherv(moved_positions.data(),
                          num_indices * NDIM,
                          MPI_DOUBLE,
                          new_positions.data(),
                          vector_counts.data(),
                          vector_offsets.data(),
                          MPI_DOUBLE,
                          IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);
    ierr = MPI_Allgatherv(moved_velocities.data(),
                          num_indices * NDIM,
                          MPI_DOUBLE,
                          new_velocities.data(),
                          vector_counts.data(),
                          vector_offsets.data(),
                          MPI_DOUBLE,
                          IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);

    // 3. Apply periodicity constraints.
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    const IntVector<NDIM> periodic_shift = grid_geom->getPeriodicShift();
    for (unsigned int k = 0; k < new_indices.size(); ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (periodic_shift[d])
            {
                double domain_length = domain_x_upper[d] - domain_x_lower[d];
                double& X = new_positions[k * NDIM + d];
                if (X < domain_x_lower[d]) X += domain_length;
                if (X >= domain_x_upper[d]) X -= domain_length;
                TBOX_ASSERT(X >= domain_x_lower[d] && X < domain_x_upper[d]);
            }
        }
    }

    // 4. Emplace the marker points.
    unsigned int num_emplaced_markers = 0;
    std::vector<bool> marker_emplaced(new_indices.size());

    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        for (MarkerPatch& marker_patch : d_marker_patches[ln])
        {
            for (unsigned int k = 0; k < new_indices.size(); ++k)
            {
                IBTK::Point position(&new_positions[k * NDIM]);
                IBTK::Vector velocity(&new_velocities[k * NDIM]);
                if (!marker_emplaced[k] && marker_patch.contains(position))
                {
                    marker_emplaced[k] = true;
                    marker_patch.insert(new_indices[k], position, velocity);
                    ++num_emplaced_markers;
                }
            }
        }
    }

    // 5. Account for points outside the computational domain which could not
    //    re-enter the domain via a periodic boundary.
    if (rank == 0)
    {
        for (unsigned int k = 0; k < new_indices.size(); ++k)
        {
            if (!marker_emplaced[k])
            {
                IBTK::Point position(&new_positions[k * NDIM]);
                // Set external velocities to zero
                IBTK::Vector velocity;
                velocity.fill(0.0);
                bool point_outside_domain = false;
                for (unsigned int d = 0; d < NDIM; ++d)
                    point_outside_domain = point_outside_domain ||
                                           ((position[d] < domain_x_lower[d]) || (domain_x_upper[d] <= position[d]));

                if (point_outside_domain)
                {
                    marker_emplaced[k] = true;
                    d_markers_outside_domain.insert(new_indices[k], position, velocity);
                    ++num_emplaced_markers;
                }
            }
        }
    }

    // 6. Check that we accounted for all points.
#ifndef NDEBUG
    ierr = MPI_Allreduce(MPI_IN_PLACE, &num_emplaced_markers, 1, MPI_UNSIGNED, MPI_SUM, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(num_emplaced_markers == new_indices.size());

    std::vector<unsigned int> marker_check(getNumberOfMarkers());
    for (int ln = d_hierarchy->getFinestLevelNumber(); ln >= 0; --ln)
    {
        for (MarkerPatch& marker_patch : d_marker_patches[ln])
        {
            for (unsigned int k = 0; k < marker_patch.size(); ++k)
            {
                const auto index = std::get<0>(marker_patch[k]);
                TBOX_ASSERT(index < long(getNumberOfMarkers()));
                marker_check[index] += 1;
            }
        }
    }
    if (rank == 0)
    {
        for (unsigned int k = 0; k < d_markers_outside_domain.size(); ++k)
        {
            const auto index = std::get<0>(d_markers_outside_domain[k]);
            TBOX_ASSERT(index < long(getNumberOfMarkers()));
            marker_check[index] += 1;
        }
    }
    ierr = MPI_Allreduce(
        MPI_IN_PLACE, marker_check.data(), marker_check.size(), MPI_UNSIGNED, MPI_SUM, IBTK_MPI::getCommunicator());
    for (unsigned int i = 0; i < getNumberOfMarkers(); ++i)
    {
        if (marker_check[i] != 1)
        {
            TBOX_ERROR(d_object_name
                       << ": Marker point " << i << " is presently owned by " << marker_check[i]
                       << " patches. The most likely cause of this error is that the CFL number is greater than 1.");
        }
    }
#endif
}
} // namespace IBTK