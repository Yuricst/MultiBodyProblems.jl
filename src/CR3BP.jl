"""Circular-Restricted Three-Body Problem"""

"""
    CR3BP(μ::Float64; DU::Float64 = 1.0, TU::Float64 = 1.0, radius_m1::Union{Nothing,Float64} = nothing, radius_m2::Union{Nothing,Float64} = nothing)

Struct for CR3BP system.
"""
mutable struct CR3BP
    μ::Float64
    μ1::Float64
    μ2::Float64
    radius_m1::Union{Nothing,Float64}
    radius_m2::Union{Nothing,Float64}
    DU::Float64
    TU::Float64
    VU::Float64

    function CR3BP(
        μ::Float64;
        DU::Float64 = 1.0,
        TU::Float64 = 1.0,
        radius_m1::Union{Nothing,Float64} = nothing,
        radius_m2::Union{Nothing,Float64} = nothing
    )
        μ1 = 1 - μ
        μ2 = μ
        VU = DU / TU
        return new(μ, μ1, μ2, radius_m1, radius_m2, DU, TU, VU)
    end
end


function Base.show(io::IO, cr3bp::CR3BP)
    println(io, "CR3BP system")
    @printf(io, "   μ         = %f\n", cr3bp.μ)
    @printf(io, "   μ1        = %f\n", cr3bp.μ1)
    @printf(io, "   μ2        = %f\n", cr3bp.μ2)
    if cr3bp.radius_m1 !== nothing
        @printf(io, "   radius_m1 = %f\n", cr3bp.radius_m1)
    end
    if cr3bp.radius_m2 !== nothing
        @printf(io, "   radius_m2 = %f\n", cr3bp.radius_m2)
    end
    @printf(io, "   DU        = %f\n", cr3bp.DU)
    @printf(io, "   TU        = %f\n", cr3bp.TU)
end


function eom_cr3bp!(dx, x, μ, t)
    # calcualte intermediate scalar quantities
    r1 = sqrt((x[1] + μ)^2 + x[2]^2 + x[3]^2)
    r2 = sqrt((x[1] - 1 + μ)^2 + x[2]^2 + x[3]^2)
    μ1_r13 = (1 - μ) / r1^3
    μ2_r23 = μ / r2^3
    # derivatives
    dx[1:3] = x[4:6]
    dx[4] =  2*x[5] + x[1] - μ1_r13*(x[1] + μ) - μ2_r23*(x[1] - 1 + μ)
    dx[5] = -2*x[4] + x[2] - μ1_r13*x[2]       - μ2_r23*x[2]
    dx[6] =                - μ1_r13*x[3]       - μ2_r23*x[3]
    return nothing
end


function dfdx_cr3bp(x, p, t)
    throw(NotImplementedError("Not implemented"))
end


function eom_stm_cr3bp!(dx_stm, x_stm, p, t)
    throw(NotImplementedError("Not implemented"))
end


function eom!(dx, x, p::CR3BP, t)
    # calcualte intermediate scalar quantities
    r1 = sqrt((x[1] + p.μ)^2 + x[2]^2 + x[3]^2)
    r2 = sqrt((x[1] - 1 + p.μ)^2 + x[2]^2 + x[3]^2)
    μ1_r13 = (1 - p.μ) / r1^3
    μ2_r23 = p.μ / r2^3
    # derivatives
    dx[1:3] = x[4:6]
    dx[4] =  2*x[5] + x[1] - μ1_r13*(x[1] + p.μ) - μ2_r23*(x[1] - 1 + p.μ)
    dx[5] = -2*x[4] + x[2] - μ1_r13*x[2]       - μ2_r23*x[2]
    dx[6] =                - μ1_r13*x[3]       - μ2_r23*x[3]
end


function eom_stm!(dx_stm, x_stm, p::CR3BP, t)
    throw(NotImplementedError("Not implemented"))
end