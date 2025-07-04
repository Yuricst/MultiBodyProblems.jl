"""Test for MultiBodyProblems.jl"""

using Test

include(joinpath(@__DIR__, "..", "src", "MultiBodyProblems.jl"))

@testset "ODEs" begin
    include("test_cr3bp.jl")
    include("test_bcr4bp.jl")
end

@testset "LPOs" begin
    include("test_jpl_lpo.jl")
    include("test_periodic_orbit.jl")
    include("test_single_shooting.jl")
end
