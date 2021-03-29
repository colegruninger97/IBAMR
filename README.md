# redbud

The MATLAB script `generate_carpet.m` can be used to set up the initial distribution of posts.
The remaining problem parameters are set in an input file (currently `input3d.NFINEST=128` is the only one in the repository).
Notice that a number of problem parameters are currently "hard coded" in this script.
If the size of the domain (`L`, `ASPECT_RATIO_X`, `ASPECT_RATIO_Y`, `ASPECT_RATIO_Z`) or the domain grid resolution (`NFINEST`) are changed in either file, they should be changed in both files.
E.g. to increase the grid resolution, update `NFINEST` in both files and regenerate the "carpet" of posts.
Notice also that `NFINEST` is a derived quantity --- it is controlled by the number of grid cells on the coarsest level of the AMR grid (`N`), the number of grid levels (`MAX_LEVELS`) and the refinement ratio between levels (`REF_RATIO`).

The post properties are set only in the MATLAB script, and can be modified independently of the rest of the problem settings (except that if the post stiffness is increased, it may be necessary to decrease the value of the maximum time step size, `DT_MAX`, in the `input3d` file).

The fluid properties are set only in the input file.
Currently they correspond to water.
