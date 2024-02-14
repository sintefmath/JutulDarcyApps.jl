using MPI, PartitionedArrays, HYPRE
using Jutul
using JutulDarcy
using LinearAlgebra
using JutulDarcyMPI


basepath = joinpath(pathof(JutulDarcy), "..", "..", "test")
for k in ["spe1", "egg", "spe9", "spe3"]
    empty!(ARGS)
    push!(ARGS, joinpath(basepath, "mrst", "$k.mat"))
    JutulDarcyMPI.julia_main()
end

empty!(ARGS)
push!(ARGS, joinpath("data", "SPE1", "BENCH_SPE1.DATA"))
JutulDarcyMPI.julia_main()

# You can add your own files here if you want certain physics / parsing to be
# included in the final result.