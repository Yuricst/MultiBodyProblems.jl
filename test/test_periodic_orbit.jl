"""Test periodic orbit"""

using LinearAlgebra
using OrdinaryDiffEq
using Test

include(joinpath(@__DIR__, "..", "src", "MultiBodyProblems.jl"))


function test_eigenspace()
    # define system
    DU = 389703.0     # km
    TU = 382981.0     # sec
    μ = 1.215058560962404e-02
    system = MultiBodyProblems.CR3BP(μ; DU=DU, TU=TU)

    # propagate initial conditions
    rv0 = [8.2812868060178335E-1, 5.4780319076101337E-28, 1.0109927300938440E-1, 1.3756588133639113E-15, 2.1658165151025582E-1, 1.9540264521098177E-15]
    period = 2.7861081325636059

    po = MultiBodyProblems.PeriodicOrbit(0.0, rv0, period)
    modes, monodromy = MultiBodyProblems.eigenspace(system, po)
    @test modes[:stable][1].value < 1.0
    @test modes[:unstable][1].value > 1.0
    @test abs(norm(modes[:oscillatory][1].value) - 1) < 1e-6
    @test abs(norm(modes[:oscillatory][2].value) - 1) < 1e-6
    @test abs(modes[:center][1].value - 1) < 1e-6
    @test abs(modes[:center][2].value - 1) < 1e-6
end


function test_invariant_manifold(;get_plot::Bool=false)
    # define system
    DU = 389703.0     # km
    TU = 382981.0     # sec
    μ = 1.215058560962404e-02
    system = MultiBodyProblems.CR3BP(μ; DU=DU, TU=TU)

    # propagate initial conditions
    rv0 = [8.2812868060178335E-1, 5.4780319076101337E-28, 1.0109927300938440E-1, 1.3756588133639113E-15, 2.1658165151025582E-1, 1.9540264521098177E-15]
    period = 2.7861081325636059

    po = MultiBodyProblems.PeriodicOrbit(0.0, rv0, period)
    modes, monodromy = MultiBodyProblems.eigenspace(system, po)
    
    t_propagate = 7.0
    sols_unstable = MultiBodyProblems.get_invariant_manifold(system, po, modes[:unstable][1], t_propagate; eps=-1e-6, N=50)
    sols_stable = MultiBodyProblems.get_invariant_manifold(system, po, modes[:stable][1], -t_propagate; eps=1e-6, N=50)
    @test true    # just check that functions ran without error

    if get_plot
        fig = Figure(size=(800,500))
        ax3d = Axis3(fig[1,1]; aspect=:data)
        for sol in sols_unstable
            lines!(ax3d, Array(sol)[1,:], Array(sol)[2,:], Array(sol)[3,:], color=:red)
        end
        for sol in sols_stable
            lines!(ax3d, Array(sol)[1,:], Array(sol)[2,:], Array(sol)[3,:], color=:blue)
        end
        display(fig)
    end
end


test_eigenspace()
test_invariant_manifold()