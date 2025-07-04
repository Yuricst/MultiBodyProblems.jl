"""Functions for generating invariant manifolds"""


"""
    get_invariant_manifold(system, po, mode, t_manifold; N=30, eps=1e-6, callback=nothing, method=Vern7(), reltol=1e-13, abstol=1e-14)

Generate an invariant manifold for a given periodic orbit and unstable mode.

# Arguments
- `system::MultiBodyProblem`: dynamics system
- `po::PeriodicOrbit`: periodic orbit 
- `mode::EigenMode`: eigenmode, must be unstable or stable
- `t_manifold::Real`: propagation time for the invariant manifold
- `N::Int`: number of initial conditions to generate
- `eps::Float64`: perturbation amplitude
- `callback`: callback function
- `method`: ODE solver method
- `reltol`: relative tolerance
- `abstol`: absolute tolerance
"""
function get_invariant_manifold(
    system::MultiBodyProblem,
    po::PeriodicOrbit,
    mode::EigenMode,
    t_manifold::Real;
    N::Int = 30,
    eps::Float64 = 1e-6,
    callback = nothing,
    method = Vern7(),
    reltol::Float64 = 1e-13,
    abstol::Float64 = 1e-14,
)
    if mode.type == :unstable
        if t_manifold <= 0.0
            warning("Given unstable mode and negative propagation time")
        end
    elseif mode.type == :stable
        if t_manifold >= 0.0
            warning("Given stable mode and positive propagation time")
        end
    else
        error("Invalid mode type $(mode.type)")
    end
    
    # propagate over one period
    prob = ODEProblem(eom_stm!, [po.x0; reshape(I(6),36)], (po.t0, po.t0 + po.period), system)
    sol = solve(prob, method, reltol=reltol, abstol=abstol)

    # evaluate starting points
    prob_branch = ODEProblem(eom!, po.x0, (0.0, 1.0), system)
    t_start = LinRange(po.t0, po.t0 + po.period, N)
    sols = ODESolution[]
    for (i,t0_branch) in enumerate(t_start)
        # get initial state & eigenvector (mapped)
        x0_lpo_i = sol(t0_branch)[1:6] 
        y = reshape(sol(t0_branch)[7:42], 6, 6)' * mode.vector
        y /= norm(y)
        
        # perturb initial state
        x0_branch = x0_lpo_i + eps * y
        
        # propagate branch
        prob_branch_i = remake(prob_branch, u0 = x0_branch, tspan = (0.0, t_manifold))
        sol_branch_i = solve(prob_branch_i, method, reltol=reltol, abstol=abstol, callback = callback)
        push!(sols, sol_branch_i)
    end
    return sols
end