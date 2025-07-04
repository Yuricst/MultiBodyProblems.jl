"""Single-shooting differential correction"""



"""
    ssdc_periodic_xzplane(p, x0, period0; kwargs...)

Single-shooting differential correction for periodic trajectory with symmetry across xz-plane
The `fix` argument specifies the degree-of-freedom to be removed from the problem.

# Arguments
- `system::MultiBodyProblem`: system struct
- `x0::Vector`: initial state
- `t0::Real`: initial time
- `period::Real`: period
- `tol_residual::Float64`: tolerance for residual
- `maxiter::Int`: maximum number of iterations
- `fix::String`: degree-of-freedom to be removed from the problem, one of: "period", "vy", "z"
- `norm_type`: type of norm to use for residual
- `verbose::Bool`: whether to print verbose output
- `method`: method to use for solving the ODE
- `reltol::Float64`: relative tolerance for the ODE solver
"""
function single_shooting_xzplane(
    system::MultiBodyProblem,
    x0::Vector,
    t0::Real,
    period::Real;
    tol_residual::Float64 = 1e-12,
    maxiter::Int = 20,
    fix::String = "period",
    norm_type = Inf,
    verbose::Bool = false,
    method = Vern7(),
    reltol::Float64 = 1e-13,
    abstol::Float64 = 1e-14,
)
    @assert fix in ["period", "vy", "z"] "fix must be one of: period, vy, z"
    @assert norm_type in [Inf, 1, 2] "norm_type must be one of: Inf, 1, 2"

    # initialize with array and period
    x0iter = deepcopy(x0)
    period_copy = deepcopy(period)
    dstatef = zeros(length(x0iter))

    # base ODE
    base_ode = ODEProblem(eom_stm!, [x0; zeros(36)], (t0, t0 + period_copy/2), system)

    # ----- iterate until convergence ----- #
    idx = 1
    converged = false
    fiters = []
    for idx in 1:maxiter        # propagate state & STM
        prob = remake(base_ode, u0=[x0iter; reshape(I(6),36)], tspan=(t0, t0 + period_copy/2))
        sol = solve(prob, method, reltol=reltol, abstol=abstol)
        statef = sol.u[end][1:length(x0iter)]
        stm = reshape(sol.u[end][length(x0iter)+1:end], (length(x0iter),length(x0iter)))';
        eom!(dstatef, statef, system, sol.t[end])
        
        # residual vector
        if length(x0iter)==4
            #        y         vx
            ferr = [ statef[2] statef[3] ]'
        elseif length(x0iter)==6
            #        y         vx        vz
            ferr = [ statef[2] statef[4] statef[6] ]'
        end

        # create free parameters vector and DF
        if fix=="vy" && length(x0iter)==4
            #      dy / dx    dvx / d(P/2)
            xi = [ x0iter[1]  period_copy/2 ]'
            df = [ [ stm[2,1] dstatef[2] ];
                   [ stm[3,1] dstatef[3] ] ];

        elseif fix=="z" && length(x0iter)==6
            #      / dx      / dvy     / d(P/2)
            xi = [ x0iter[1] x0iter[5] period_copy/2 ]'
            df = [ [ stm[2,1] stm[2,5] dstatef[2] ];
                   [ stm[4,1] stm[4,5] dstatef[4] ];
                   [ stm[6,1] stm[6,5] dstatef[6] ] ];

        elseif fix=="period" && length(x0iter)==4
            #      dy / dx    dvx / dvy
            xi = [ x0iter[1]  x0iter[4] ]'
            df = [ [ stm[2,1] stm[2,4] ];
                   [ stm[3,1] stm[3,4] ] ];

        elseif fix=="period" && length(x0iter)==6
            #      / dx      / dz      / dvy
            xi = [ x0iter[1] x0iter[3] x0iter[5] ]'
            df = [ [ stm[2,1] stm[2,3] stm[2,5] ];
                   [ stm[4,1] stm[4,3] stm[4,5] ];
                   [ stm[6,1] stm[6,3] stm[6,5] ] ];
        end

        # correct state
        xii = xi - inv(df)*ferr
        push!(fiters, norm(ferr))

        # update state
        if fix=="vy" && length(x0iter)==4
            #      x         P/2
            x0iter[1] = xii[1];
            period_copy = 2xii[2];
        elseif fix=="z" && length(x0iter)==6
            #      x         vy        P/2
            x0iter[1] = xii[1];
            x0iter[5] = xii[2];
            period_copy = 2xii[3];
        elseif fix=="period" && length(x0iter)==4
            #      x         vy
            x0iter[1] = xii[1];
            x0iter[4] = xii[2];
        elseif fix=="period" && length(x0iter)==6
            #      x         z         vy
            x0iter[1] = xii[1];
            x0iter[3] = xii[2];
            x0iter[5] = xii[3];
        end

        # check convergence
        if verbose
            @printf("Iteration %i: residual: %s\n", idx, norm(ferr, norm_type))
        end
        if norm(ferr, norm_type) <= tol_residual
            if verbose
                @printf("Converged at iteration %i\n", idx)
            end
            converged = true
            break
        end
    end
    if converged
        po = PeriodicOrbit(t0, x0iter, period_copy)
    else
        po = nothing
    end
    return po, converged, fiters
end
