"""Test for CR3BP system"""

using LinearAlgebra
using OrdinaryDiffEq
using Test

include(joinpath(@__DIR__, "..", "src", "MultiBodyProblems.jl"))


function test_cr3bp_eom()
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


function test_cr3bp_dfdx()
    # define system
    DU = 389703.0     # km
    μ = 1.215058560962404e-02
    system = MultiBodyProblems.CR3BP(μ; DU=DU)

    # define initial conditions
    rv0 = [1.0809931218390707, 0.03, -2.0235953267405354E-01, -0.03, -1.9895001215078018E-01, 0.02]

    # numerically compute dfdx
    A_num = zeros(6,6)
    h = 1e-8
    for i in 1:6
        x_plus = copy(rv0)
        x_plus[i] += h
        x_minus = copy(rv0)
        x_minus[i] -= h
        dx_plus = zeros(6)
        dx_minus = zeros(6)
        MultiBodyProblems.eom!(dx_plus, x_plus, system, 0.0)
        MultiBodyProblems.eom!(dx_minus, x_minus, system, 0.0)
        A_num[:,i] = (dx_plus - dx_minus) / (2*h)
    end

    # compare with analytical dfdx
    A = MultiBodyProblems.dfdx(rv0, system, 0.0)
    @test norm(A_num - A) < 1e-7

    A_explicit = MultiBodyProblems.dfdx_cr3bp(rv0, system.μ, 0.0)
    @test norm(A_num - A_explicit) < 1e-7
end


function test_cr3bp_stm()
    # define system
    DU = 389703.0     # km
    μ = 1.215058560962404e-02
    system = MultiBodyProblems.CR3BP(μ; DU=DU)

    # define initial conditions
    rv0 = [1.0809931218390707, 0.0, -2.0235953267405354E-01, 0.0, -1.9895001215078018E-01, 0.0]
    period = 2.3538670417546639

    # solve ODE over one period
    x0_stm = [rv0; reshape(I(6), 36)]
    prob = ODEProblem(MultiBodyProblems.eom_stm!, x0_stm, (0.0, period), system)   # eom implied by parameters
    sol = solve(prob, Vern7(), reltol=1e-13, abstol=1e-14)
    STM_1 = reshape(sol.u[end][7:42], 6, 6)'
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test sol.u[end][1:6] ≈ rv0 atol=1e-12

    prob_2 = ODEProblem(MultiBodyProblems.eom_stm_cr3bp!, x0_stm, (0.0, period), system.μ)   # eom_cr3bp! requires μ as a parameter
    sol_2 = solve(prob_2, Vern7(), reltol=1e-13, abstol=1e-14)
    STM_2 = reshape(sol_2.u[end][7:42], 6, 6)'
    @test sol_2.retcode == SciMLBase.ReturnCode.Success
    @test sol_2.u[end][1:6] ≈ rv0 atol=1e-12

    @test STM_1 == STM_2
    @test sol.u[end] == sol_2.u[end]

    # numerically compute STM
    STM_num = zeros(6,6)
    h = 1e-8
    for i in 1:6
        x0_plus = copy(rv0)
        x0_plus[i] += h
        x0_minus = copy(rv0)
        x0_minus[i] -= h
        prob_plus = ODEProblem(MultiBodyProblems.eom!, x0_plus, (0.0, period), system)
        prob_minus = ODEProblem(MultiBodyProblems.eom!, x0_minus, (0.0, period), system)
        sol_plus = solve(prob_plus, Vern7(), reltol=1e-13, abstol=1e-14)
        sol_minus = solve(prob_minus, Vern7(), reltol=1e-13, abstol=1e-14)
        STM_num[:,i] = (sol_plus.u[end] - sol_minus.u[end]) / (2*h)
    end
    @test maximum(abs.(STM_1 - STM_num)) < 1e-6
    @test maximum(abs.(STM_2 - STM_num)) < 1e-6
end


test_cr3bp_eom()
test_cr3bp_dfdx()
test_cr3bp_stm()