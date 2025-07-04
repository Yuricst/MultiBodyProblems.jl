"""Test for JPL LPOs"""

using LinearAlgebra
using OrdinaryDiffEq
using Test

include(joinpath(@__DIR__, "..", "src", "MultiBodyProblems.jl"))


function test_jpl_lpo()
    loadfile = joinpath(@__DIR__, "..", "data", "jpl_lpos", "EarthMoon_halo_L2S.json")
    dict_family, df = MultiBodyProblems.read_jpl_threebody_periodic_orbits(loadfile)
    @test true
end

test_jpl_lpo()