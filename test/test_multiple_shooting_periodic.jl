"""Test single shooting"""

using LinearAlgebra
using OrdinaryDiffEq
using Test

include(joinpath(@__DIR__, "..", "src", "MultiBodyProblems.jl"))


function test_multiple_shooting_periodic(;get_plot::Bool = false)
    # define system
    DU = 389703.0     # km
    TU = 382981.0     # sec
    μ = 1.215058560962404e-02
    system = MultiBodyProblems.CR3BP(μ; DU=DU, TU=TU)

    # propagate initial conditions
    rv0 = [1.1694440291316723, 0.0, 9.6468484190090548E-2, 0.0, -1.9375504587745879E-1, 0.0]
    period = 3.3   # 3.3341632624577340

    # propagate initial guess
    prob_ig = ODEProblem(MultiBodyProblems.eom!, rv0, (-period/2, period/2), system)   # eom implied by parameters
    sol_ig = solve(prob_ig, Vern8(), reltol=1e-14, abstol=1e-14)
    times = Vector(LinRange(-period/2, period/2, 7))
    xs = hcat(sol_ig(times[1:end-1])...)

    po, converged, DF, sols = MultiBodyProblems.multiple_shooting_periodic(
        system, xs, times; tol_residual=1e-11, alpha=0.9, maxiter=200, verbose=true
    )
    @test converged == true

    # propagate 
    prob = ODEProblem(MultiBodyProblems.eom!, po.x0, (-period/2, period/2), system)   # eom implied by parameters
    sol = solve(prob, Vern8(), reltol=1e-14, abstol=1e-14)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    #@show maximum(abs.(sol.u[end] - po.x0))
    @test maximum(abs.(sol.u[end] - po.x0)) < 1e-10   # via single propagation

    # plot
    if get_plot
        fig = Figure(size=(500,500))
        ax3d = Axis3(fig[1,1]; aspect=:data)
        for sol in sols
            lines!(ax3d, Array(sol)[1,:], Array(sol)[2,:], Array(sol)[3,:])
        end
        lines!(ax3d, Array(sol_ig)[1,:], Array(sol_ig)[2,:], Array(sol_ig)[3,:], color=:red)
        display(fig)
    end
    return sparse(DF)
end


DF = test_multiple_shooting_periodic(;get_plot = false)