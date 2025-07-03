"""Periodic orbit struct"""


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