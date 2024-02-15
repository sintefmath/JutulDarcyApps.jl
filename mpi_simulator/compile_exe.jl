using PackageCompiler, Jutul, JutulDarcy
# This can take a while.
create_app("JutulDarcyMPI", "compiled_simulator",
    precompile_execution_file = "precompile_jutul_darcy_mpi.jl", # Precompilation script
    force = true,                                                # Delete existing files
    incremental=true,                                            # Add onto existing Julia sysimage
    sysimage_build_args = Cmd(["-O2"])                           # Set Julia flags for precompilation
)
