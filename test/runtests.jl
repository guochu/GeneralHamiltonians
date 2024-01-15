push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/DMRG/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/InfiniteDMRG/src")
push!(LOAD_PATH, dirname(Base.@__DIR__) * "/src")

using Test, Random
using SphericalTensors, DMRG, InfiniteDMRG
using GeneralHamiltonians

Random.seed!(1342)

include("util.jl")

# include("auxiliary.jl")


# include("finitesystem.jl")

include("infinitesystem.jl")
