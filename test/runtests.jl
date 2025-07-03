"""Test for MultiBodyProblems.jl"""

using Test

include(joinpath(@__DIR__, "..", "src", "MultiBodyProblems.jl"))

@testset "ODEs" begin
    include("test_cr3bp.jl")
end
