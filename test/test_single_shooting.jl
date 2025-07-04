"""Test single shooting"""

using LinearAlgebra
using OrdinaryDiffEq
using Test

include(joinpath(@__DIR__, "..", "src", "MultiBodyProblems.jl"))


function test_single_shooting(;get_plot::Bool = false)
    # define system
    DU = 389703.0     # km
    TU = 382981.0     # sec
    μ = 1.215058560962404e-02
    system = MultiBodyProblems.CR3BP(μ; DU=DU, TU=TU)

    # propagate initial conditions
    rv0 = [8.2812868060178335E-1, 5.4780319076101337E-28, 1.0109927300938440E-1, 1.3756588133639113E-15, 2.1658165151025582E-1, 1.9540264521098177E-15]
    period = 2.76

    po, converged, fiters = MultiBodyProblems.single_shooting_xzplane(system, rv0, 0.0, period; verbose=false)
    @test converged == true

    # propagate 
    prob = ODEProblem(MultiBodyProblems.eom!, po.x0, (0.0, period), system)   # eom implied by parameters
    sol = solve(prob, Vern7(), reltol=1e-13, abstol=1e-14)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test sol.u[end] ≈ po.x0 atol=1e-12

    # plot
    if get_plot
        fig = Figure(size=(800,500))
        ax3d = Axis3(fig[1,1]; aspect=:data)
        lines!(ax3d, Array(sol)[1,:], Array(sol)[2,:], Array(sol)[3,:], color=:red)
        display(fig)
    end
end


test_single_shooting(;get_plot = false)