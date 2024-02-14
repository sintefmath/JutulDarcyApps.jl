# Building a stand-alone JutulDarcy.jl reservoir simulator

This folder contains the necessary files to build a standalone reservoir simulator from [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl) using [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl).

The main motivation is to create applications for users and systems where Julia is not installed, integration in workflows where compilation is impractical and for users who prefer to run simulations from the command line.

## Building the application

### Setting correct directory and project

```bash
cd mpi_simulator/
julia --project=.
```

### First time setup

If it is your first time building the simulator, execute the following to add the local `JutulDarcyMPI` package to the environment:

```julia
] # Enter package mode
dev ./JutulDarcyMPI # develop local package that contains application entry point
```

### Finally: Build the application

If this is your first time building the application, this next step will take some time: JutulDarcy and all other application dependencies are downloaded, precompiled and finally compiled. During the compilation process, a few standard test cases are run.

```julia
include("compile_exe.jl")
```

The output should end with something that looks similar to the following:

```bash
[ Info: PackageCompiler: Done
✔ [10m:42s] PackageCompiler: compiling incremental system image
```

If PackageCompiler reports success, the application has been successfully built for your platform.

## Usage of application

The built application will appear in the subdirectory `mpi_simulator/compiled_similator`. The entire `compiled_simulator` folder and all subdirectories (`bin`, `lib`, etc) of that folder are required to run the simulator. This is important if you plan on moving the executable to another machine.

### Running a case in serial

You can now try to run the included test case:

#### Windows style paths

```bash
.\compiled_simulator\bin\JutulDarcyMPI.exe .\data\SPE1\BENCH_SPE1.DATA
```

#### Linux and Mac OS style paths

```bash
./compiled_simulator/bin/JutulDarcyMPI ./data/SPE1/BENCH_SPE1.DATA
```

The case should run in less than a second on most computers, indicating that the simulator does not need to compile anything.

### Running .DATA files

You can run a subset of .DATA files. This is done using the [GeoEnergyIO.jl parser](https://github.com/sintefmath/GeoEnergyIO.jl). Not all keywords are supported. If you have an example you would like to get working, post it at the issue page there.

### Running MRST (.MAT) files

The compiled simulator can also read cases exported from [MRST](https://www.sintef.no/projectweb/mrst/) using the [jutul module](https://github.com/SINTEF-AppliedCompSci/MRST/tree/main/model-io/jutul). There is no special option to run these files as the format is automatically detected from the file extension.

### Reading the results

The results are currently written as [JLD2](https://github.com/JuliaIO/JLD2.jl) files. This format is a subset of HDF5 and can be read using standard HDF5 readers. You can recover the simulation results in Julia session using the `read_results` function:

```julia
using Jutul
states, reports = read_results("/path/to/output/folder")
```

Adding more user friendly output (e.g. well solutions as csv files) is a future goal.

### Command line options

There are many options to that can be used to customize tolerances, linear solvers, time-stepping and output level:

```bash
.\compiled_simulator\bin\JutulDarcyMPI.exe
required argument filename was not provided
usage: <PROGRAM> [--tol-cnv TOL_CNV] [--tol-mb TOL_MB]
                 [--tol-dp-well TOL_DP_WELL]
                 [--multisegment-wells MULTISEGMENT_WELLS]
                 [--relaxation RELAXATION]
                 [--tol-factor-final-iteration TOL_FACTOR_FINAL_ITERATION]
                 [--max-nonlinear-iterations MAX_NONLINEAR_ITERATIONS]
                 [--max-timestep-cuts MAX_TIMESTEP_CUTS] [--rtol RTOL]
                 [--backend BACKEND] [--linear-solver LINEAR_SOLVER]
                 [--precond PRECOND] [--timesteps TIMESTEPS]
                 [--timesteps-target-ds TIMESTEPS_TARGET_DS]
                 [--timesteps-target-iterations TIMESTEPS_TARGET_ITERATIONS]
                 [--max-timestep MAX_TIMESTEP]
                 [--initial-timestep INITIAL_TIMESTEP]
                 [--info-level INFO_LEVEL] [--verbose VERBOSE]
                 [--output-path OUTPUT_PATH] filename
```

You can get more information on the different options by passing `--help`.

### Running using MPI / parallel computing

The executable is built with MPI support for parallel execution. To run using MPI, you must call the function using the correct MPI executable. For more information on this topic, see [the Julia MPI documentation](https://juliaparallel.org/MPI.jl/stable/). If the wrong MPI executable is used, each MPI process will simulate the case separately - leading to slowdown instead of speedup as redundant work is performed.

The easiest way to get started is to install the default Julia MPI executable using the MPI module that is already installed in your build environment:

```julia
julia> using MPI
julia> MPI.install_mpiexecjl()
julia> MPI.mpiexec() # Will display the path
```

You will then get a listing that shows you the path of the MPI executable. In my case, that path is `C:\Users\olavm\.julia\artifacts\4a195fc1f9aa324d571de9f2e6efdf6bb3562edf\bin\mpiexec.exe`. I can then use that executable to simulate a case in parallel.

For example, to simulate with 2 processes I could do the following:

```bash
C:\Users\olavm\.julia\artifacts\4a195fc1f9aa324d571de9f2e6efdf6bb3562edf\bin\mpiexec.exe -n 2 .\compiled_simulator\bin\JutulDarcyMPI.exe data\SPE1\BENCH_SPE1.DATA
```

The SPE1 model is tiny and used to demonstrate the concept. Do not expect any performance benefit from using parallel computing for this case.

## Advanced usage

The `precompile_jutul_darcy_mpi.jl` script is used to cover code paths for compilation. If your compiled simulator takes a bit too long to start up, you can add additional cases in this file. These cases are run during the compilation stage and allows you to cover more code paths. The downside is that the cases will then be run during compilation, which may slow things down if you add large or complex cases.

## Contributing

This application was primarily written as a demonstrator for compilation workflows with JutulDarcy. Contributions to ease of use, functionality and documentation are welcome.
