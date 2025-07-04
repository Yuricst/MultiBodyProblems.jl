"""Multiple-shooting differential correction"""


function multiple_shooting_periodic(
    system::MultiBodyProblem,
    xs,
    times::Vector;
    tol_residual::Float64 = 1e-11,
    alpha::Float64 = 0.75,
    maxiter::Int = 200,
    norm_type = Inf,
    verbose::Bool = true,
    method = Vern8(),
    reltol::Float64 = 1e-14,
    abstol::Float64 = 1e-14,
)
    N = length(times)
    Nseg = N - 1
    @assert size(xs) == (6,Nseg) "xs must be a matrix of size (6,Nseg)"

    # initialize storages & ODEProblem
    prob_base = ODEProblem(eom_stm!, [xs[:,1]; zeros(36)], (0.0,0.1), system)
    F = zeros(6*Nseg)
    xs_flat = reshape(xs,6*Nseg)
    converged = false
    sols = [solve(prob_base, method, reltol=reltol, abstol=abstol) for _ in 1:Nseg]  # prepare storage

    # initialize Jacobian
    DF = zeros(6*Nseg, 6*Nseg)
    for j_seg in 1:Nseg
        if j_seg < Nseg
            DF[1+6*(j_seg-1):6*j_seg, 1+6*j_seg:6*(j_seg+1)] = I(6)
        else
            DF[1+6*(j_seg-1):6*j_seg, 1:6] = I(6)
        end
    end

    for idx in 1:maxiter
        # propagate each segment
        for j_seg in 1:Nseg
            prob = remake(prob_base, u0=[xs_flat[1+6*(j_seg-1):6*j_seg]; reshape(I(6),36)], tspan=(times[j_seg], times[j_seg+1]))
            sol = solve(prob, method, reltol=reltol, abstol=abstol)
            
            DF[1+6*(j_seg-1):6*j_seg, 1+6*(j_seg-1):6*j_seg] = -reshape(sol.u[end][7:42],6,6)'
            if j_seg < Nseg
                F[1+6*(j_seg-1):6*j_seg] = xs_flat[1+6*j_seg:6*(j_seg+1)] - sol.u[end][1:6]
            else
                F[1+6*(j_seg-1):6*j_seg] = xs_flat[1:6] - sol.u[end][1:6]
            end

            sols[j_seg] = sol
        end

        # check convergence
        Ferr = norm(F, norm_type)
        if verbose
            @printf("Iteration %i: residual: %1.4e\n", idx, Ferr)
        end
        if Ferr <= tol_residual
            if verbose
                @printf("Converged at iteration %i\n", idx)
            end
            converged = true
            break
        end

        # compute update
        xs_flat -= alpha * inv(DF)*F
    end
    if converged
        po = PeriodicOrbit(times[1], xs_flat[1:6], times[end] - times[1])
    else
        po = nothing
    end
    return po, converged, DF, sols
end