"""Periodic orbit struct"""


mutable struct PeriodicOrbit
    t0::Real
    x0::Vector
    period::Real

    function PeriodicOrbit(t0::Real, x0::Vector, period::Real)
        return new(t0, x0, period)
    end
end
