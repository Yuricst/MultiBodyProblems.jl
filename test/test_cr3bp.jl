"""Test for CR3BP system"""

using OrdinaryDiffEq
using Test

include(joinpath(@__DIR__, "..", "src", "MultiBodyProblems.jl"))


function test_cr3bp()
    # define system
    DU = 389703.0     # km
    TU = 382981.0     # sec
    μ = 1.215058560962404e-02
    system = MultiBodyProblems.CR3BP(μ; DU=DU, TU=TU)

    # propagate initial conditions
    rv0 = [1.0809931218390707, 0.0, -2.0235953267405354E-01, 0.0, -1.9895001215078018E-01, 0.0]
    period = 2.3538670417546639

    # solve ODE over one period
    prob = ODEProblem(MultiBodyProblems.eom!, rv0, (0.0, period), system)   # eom implied by parameters
    sol = solve(prob, Vern7(), reltol=1e-13, abstol=1e-14)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test sol.u[end] ≈ rv0 atol=1e-12

    prob_2 = ODEProblem(MultiBodyProblems.eom_cr3bp!, rv0, (0.0, period), system.μ)   # eom_cr3bp! requires μ as a parameter
    sol_2 = solve(prob_2, Vern7(), reltol=1e-13, abstol=1e-14)
    @test sol_2.retcode == SciMLBase.ReturnCode.Success
    @test sol_2.u[end] ≈ rv0 atol=1e-12
end

test_cr3bp()