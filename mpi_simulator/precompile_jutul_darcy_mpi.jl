using MPI, PartitionedArrays, HYPRE
using Jutul
using JutulDarcy
using LinearAlgebra
using JutulDarcyMPI

pth, = splitdir(pathof(JutulDarcy))
basepath = joinpath(pth, "..", "test")
for k in ["spe1", "egg", "spe9", "spe3"]
    empty!(ARGS)
    push!(ARGS, joinpath(basepath, "mrst", "$k.mat"))
    if isfile(ARGS[1])
        JutulDarcyMPI.julia_main()
    end
end

empty!(ARGS)
push!(ARGS, joinpath("data", "SPE1", "BENCH_SPE1.DATA"))
if isfile(ARGS[1])
    JutulDarcyMPI.julia_main()
end

# You can add your own files here if you want certain physics / parsing to be
# included in the final result.
