using MPI, PartitionedArrays, HYPRE
using Jutul
using JutulDarcy
using LinearAlgebra
using JutulDarcyMPI

compile_nldd = true
##
pth, = splitdir(pathof(JutulDarcy))
basepath = joinpath(pth, "..", "test")
for k in ["spe1", "egg", "spe9", "spe3"]
    empty!(ARGS)
    push!(ARGS, joinpath(basepath, "mrst", "$k.mat"))
    push!(ARGS, "--number-of-steps=1")
    if isfile(ARGS[1])
        JutulDarcyMPI.julia_main()
        if compile_nldd
            push!(ARGS, "--method=nldd")
            push!(ARGS, "--nldd-solve-tol-mobility=0.01")
            push!(ARGS, "--nldd-solve-tol-composition=0.2")
            push!(ARGS, "--nldd-solve-tol-saturations=0.2")
            JutulDarcyMPI.julia_main()
        end
    end
end
##
empty!(ARGS)
push!(ARGS, joinpath("data", "SPE1", "BENCH_SPE1.DATA"))
if isfile(ARGS[1])
    # Compile verbosity etc
    push!(ARGS, "--number-of-steps=1")
    push!(ARGS, "--verbose")
    push!(ARGS, "--info-level=4")
    JutulDarcyMPI.julia_main()
end

# You can add your own files here if you want certain physics / parsing to be
# included in the final result.
