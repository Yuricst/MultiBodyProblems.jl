"""Test for CR3BP system"""

using LinearAlgebra
using OrdinaryDiffEq
using Test

include(joinpath(@__DIR__, "..", "src", "MultiBodyProblems.jl"))


function test_bcr4bp_eom()
    # define system
    DU = 389703.0     # km
    TU = 382981.0     # sec
    μ = 0.012150584269940354
    μ_s = 328900.5598102475
    a_s = 389.1779396462019
    ω_s = -0.9253018261815922
    system = MultiBodyProblems.BCR4BP(μ, μ_s, a_s, ω_s; DU=DU, TU=TU)

    # propagate initial conditions
    x0 = [0.9853006895390061, 0.0021723003109550003, 0.0034046027917,
      -1.198746524193054, 1.038423839870351, -1.559959799100589]
    
    tspan = (0.0, -22.970560065381424)

    # solve ODE over one period
    prob = ODEProblem(MultiBodyProblems.eom!, x0, tspan, system)   # eom implied by parameters
    sol = solve(prob, Vern7(), reltol=1e-13, abstol=1e-14)
    @test sol.retcode == SciMLBase.ReturnCode.Success

    prob_2 = ODEProblem(MultiBodyProblems.eom_bcr4bp!, x0, tspan, [system.μ, system.μ_s, system.a_s, system.ω_s])   # eom_cr3bp! requires μ as a parameter
    sol_2 = solve(prob_2, Vern7(), reltol=1e-13, abstol=1e-14)
    @test sol_2.retcode == SciMLBase.ReturnCode.Success
    
    xf_test = [0.3810885520047902, -0.49442795718815136, 0.046195317004687514,
              -1.445010410675033, 0.8597808943242603, 0.03253211976324759]
    @test maximum(abs.(xf_test - sol.u[end][1:6])) < 1e-9     # low tolerance needed for github actions
    @test maximum(abs.(xf_test - sol_2.u[end][1:6])) < 1e-9   # low tolerance needed for github actions

    # # figure
    # fig = Figure(size=(600,600))
    # # ax3d = Axis3(fig[1,1], aspect=:data)
    # # lines!(ax3d, Array(sol)[1,:], Array(sol)[2,:], Array(sol)[3,:])
    # # lines!(ax3d, Array(sol_2)[1,:], Array(sol_2)[2,:], Array(sol_2)[3,:])
    # ax = Axis(fig[1,1],)
    # lines!(ax, Array(sol)[1,:], Array(sol)[2,:])
    # display(fig)
end


function test_bcr4bp_dfdx()
    # define system
    DU = 389703.0     # km
    TU = 382981.0     # sec
    μ = 0.012150584269940354
    μ_s = 328900.5598102475
    a_s = 389.1779396462019
    ω_s = -0.9253018261815922
    system = MultiBodyProblems.BCR4BP(μ, μ_s, a_s, ω_s; DU=DU, TU=TU)

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

    A_explicit = MultiBodyProblems.dfdx_bcr4bp(rv0, [system.μ, system.μ_s, system.a_s, system.ω_s], 0.0)
    @test norm(A_num - A_explicit) < 1e-7
end


function test_bcr4bp_stm()
    # define system
    DU = 389703.0     # km
    TU = 382981.0     # sec
    μ = 0.012150584269940354
    μ_s = 328900.5598102475
    a_s = 389.1779396462019
    ω_s = -0.9253018261815922
    system = MultiBodyProblems.BCR4BP(μ, μ_s, a_s, ω_s; DU=DU, TU=TU)

    # propagate initial conditions
    x0 = [0.9853006895390061, 0.0021723003109550003, 0.0034046027917,
      -1.198746524193054, 1.038423839870351, -1.559959799100589]
    
    tspan = (0.0, -0.5)

    # # solve ODE over one period
    x0_stm = [x0; reshape(I(6), 36)]
    prob = ODEProblem(MultiBodyProblems.eom_stm!, x0_stm, tspan, system)   # eom implied by parameters
    sol = solve(prob, Vern7(), reltol=1e-13, abstol=1e-14)
    STM_1 = reshape(sol.u[end][7:42], 6, 6)'
    @test sol.retcode == SciMLBase.ReturnCode.Success

    prob_2 = ODEProblem(MultiBodyProblems.eom_stm_bcr4bp!, x0_stm, tspan, [system.μ, system.μ_s, system.a_s, system.ω_s])
    sol_2 = solve(prob_2, Vern7(), reltol=1e-13, abstol=1e-14)
    STM_2 = reshape(sol_2.u[end][7:42], 6, 6)'
    @test sol_2.retcode == SciMLBase.ReturnCode.Success

    @test STM_1 == STM_2
    @test sol.u[end] == sol_2.u[end]

    # numerically compute STM
    STM_num = zeros(6,6)
    h = 1e-8
    for i in 1:6
        x0_plus = copy(x0)
        x0_plus[i] += h
        x0_minus = copy(x0)
        x0_minus[i] -= h
        prob_plus = ODEProblem(MultiBodyProblems.eom!, x0_plus, tspan, system)
        prob_minus = ODEProblem(MultiBodyProblems.eom!, x0_minus, tspan, system)
        sol_plus = solve(prob_plus, Vern7(), reltol=1e-13, abstol=1e-14)
        sol_minus = solve(prob_minus, Vern7(), reltol=1e-13, abstol=1e-14)
        STM_num[:,i] = (sol_plus.u[end] - sol_minus.u[end]) / (2*h)
    end
    @test maximum(abs.(STM_1 - STM_num)) < 3e-5
    @test maximum(abs.(STM_2 - STM_num)) < 3e-5
end

test_bcr4bp_eom()
test_bcr4bp_dfdx()
test_bcr4bp_stm()