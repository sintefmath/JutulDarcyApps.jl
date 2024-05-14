using PackageCompiler, Jutul, JutulDarcy
# This can take a while.
cd(@__DIR__)
@time create_app("JutulDarcyMPI", "compiled_simulator",
    precompile_execution_file = "precompile_jutul_darcy_mpi.jl", # Precompilation script
    force = true,                                                # Delete existing files
    incremental=true,                                            # Add onto existing Julia sysimage
    sysimage_build_args = Cmd(["-O2"])                           # Set Julia flags for precompilation
)
##
using MPI
mpi_exec_path = only(MPI.mpiexec().exec)
if Sys.iswindows()
    exe_name = "JutulDarcyMPI.exe"
else
    exe_name = "JutulDarcyMPI"
end
exe_path = joinpath(@__DIR__, "compiled_simulator", "bin", exe_name)
println("MPI simulator compiled\n\t$exe_path")
println("MPI executable path:\n\t$mpi_exec_path")

test_file_path = joinpath(@__DIR__, "data", "SPE1", "BENCH_SPE1.DATA")
println("\nTo run simulation in serial:\n")
println("$exe_path $test_file_path")

println("\nTo run simulation in parallel as e.g. two processes:\n")
println("$mpi_exec_path -n 2 $exe_path $test_file_path")
##
println("Listing help page:")
run(`$exe_path --help`)
