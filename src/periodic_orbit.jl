"""Periodic orbit struct"""


"""
    dict_family, df = read_jpl_threebody_periodic_orbits(loadfile::String)

Load periodic orbits data from json file from JPL database

# Returns
- `dict_family::Dict`: dictionary containing the family name and other metadata
- `df::DataFrame`: DataFrame containing the periodic orbits data
"""
function read_jpl_threebody_periodic_orbits(loadfile::String)
    # get initial states and periods from file
    dict_family = JSON.parsefile(loadfile)
    df = DataFrame(parse.(Float64,hcat(dict_family["data"]...))', dict_family["fields"])
    return dict_family, df
end


"""
    PeriodicOrbit(t0::Real, x0::Vector, period::Real; family_name::String = "Unspecified")

Periodic orbit struct

# Arguments
- `t0::Real`: initial time
- `x0::Vector`: initial state
- `period::Real`: period
"""
mutable struct PeriodicOrbit
    family_name::String
    t0::Real
    x0::Vector
    period::Real

    function PeriodicOrbit(t0::Real, x0::Vector, period::Real; family_name::String = "Unspecified")
        return new(family_name, t0, x0, period)
    end
end


function Base.show(io::IO, po::PeriodicOrbit)
    println("PeriodicOrbit struct")
    @printf("  family : %s\n", po.family_name)
    @printf("  t0     : %f\n", po.t0)
    @printf("  period : %f\n", po.period)
    @printf("  r0     : % 1.4e, % 1.4e, % 1.4e\n", po.x0[1:3]...)
    @printf("  v0     : % 1.4e, % 1.4e, % 1.4e\n", po.x0[4:6]...)
end


function eval_periodic_solution(sol::ODESolution, t::Real, period::Real)
    return sol(mod(t, period))
end


mutable struct EigenMode
    type::Symbol
    value
    vector::Vector
end


function Base.show(io::IO, em::EigenMode)
    println("EigenMode struct")
    @printf("  type  : %s\n", em.type)
    @printf("  value : % 1.4e + % 1.4ei\n", real(em.value), imag(em.value))
end



function eigenspace(
    system::MultiBodyProblem,
    po::PeriodicOrbit;
    method = Vern7(),
    reltol::Float64 = 1e-13,
    abstol::Float64 = 1e-14,
    eps = 1e-6,
)
    # propagate over one period
    prob = ODEProblem(eom_stm!, [po.x0; reshape(I(6),36)], (po.t0, po.t0 + po.period), system)
    sol = solve(prob, method, reltol=reltol, abstol=abstol)

    # compute monodromy matrix
    STM = reshape(sol.u[end][7:42], 6, 6)'

    # compute eigenvalues and eigenvectors
    eig_res = eigen(STM)
    eigvals = eig_res.values
    eigvecs = eig_res.vectors

    # construct modes
    modes = Dict(
        :stable => [],
        :unstable => [],
        :oscillatory => [],
        :center => [],
    )
    for (i,λ) in enumerate(eigvals)
        if abs(imag(λ)) < eps
            if abs(real(λ) - 1.0) < eps
                push!(modes[:center], EigenMode(:center, real(λ), real(eigvecs[:,i])))
            elseif real(λ) < 1.0
                push!(modes[:stable], EigenMode(:stable, real(λ), real(eigvecs[:,i])))
            else
                push!(modes[:unstable], EigenMode(:unstable, real(λ), real(eigvecs[:,i])))
            end
        else
            push!(modes[:oscillatory], EigenMode(:oscillatory, λ, eigvecs[:,i]))
        end
    end
    return modes, STM
end