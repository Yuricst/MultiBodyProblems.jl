"""Circular-Restricted Three-Body Problem"""

"""
    CR3BP(μ::Float64; DU::Float64 = 1.0, TU::Float64 = 1.0, radius_m1::Union{Nothing,Float64} = nothing, radius_m2::Union{Nothing,Float64} = nothing)

Struct for CR3BP system.
"""
mutable struct CR3BP <: MultiBodyProblem
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


"""
    eom_cr3bp!(dx, x, μ, t)

EOM for CR3BP system compatible with `DifferentialEquations.jl` in mutating form.
"""
function eom_cr3bp!(dx, x, μ, t)
    # calcualte intermediate scalar quantities
    r1 = sqrt((x[1] + μ)^2 + x[2]^2 + x[3]^2)
    r2 = sqrt((x[1] - 1 + μ)^2 + x[2]^2 + x[3]^2)
    μ1_r13 = (1 - μ) / r1^3
    μ2_r23 = μ / r2^3
    # derivatives of state
    dx[1:3] = x[4:6]
    dx[4] =  2*x[5] + x[1] - μ1_r13*(x[1] + μ) - μ2_r23*(x[1] - 1 + μ)
    dx[5] = -2*x[4] + x[2] - μ1_r13*x[2]       - μ2_r23*x[2]
    dx[6] =                - μ1_r13*x[3]       - μ2_r23*x[3]
    return nothing
end


"""
    dfdx_cr3bp(x, μ, t)

Compute Jacobian of CR3BP EOM.
"""
function dfdx_cr3bp(x, μ, t)
    A = zeros(6,6)
    A[1:3,4:6] = I(3)
    A[4,5] = 2.0
    A[5,4] = -2.0
    G1 = gravity_gradient_matrix(1-μ, [x[1] + μ, x[2], x[3]])
    G2 = gravity_gradient_matrix(μ, [x[1] - 1 + μ, x[2], x[3]])
    A[4:6,1:3] = [1 0 0; 0 1 0; 0 0 0] + G1 + G2
    return A
end


"""
    eom_stm_cr3bp!(dx_stm, x_stm, μ, t)

EOM for CR3BP system with STM compatible with `DifferentialEquations.jl` in mutating form.
The state `x_stm` consists of the Cartesian r, v vectors & STM flattened in row-major order.
"""
function eom_stm_cr3bp!(dx_stm, x_stm, μ, t)
    # calcualte intermediate scalar quantities
    r1 = sqrt((x_stm[1] + μ)^2 + x_stm[2]^2 + x_stm[3]^2)
    r2 = sqrt((x_stm[1] - 1 + μ)^2 + x_stm[2]^2 + x_stm[3]^2)
    μ1_r13 = (1 - μ) / r1^3
    μ2_r23 = μ / r2^3
    # derivatives of state
    dx_stm[1:3] = x_stm[4:6]
    dx_stm[4] =  2*x_stm[5] + x_stm[1] - μ1_r13*(x_stm[1] + μ) - μ2_r23*(x_stm[1] - 1 + μ)
    dx_stm[5] = -2*x_stm[4] + x_stm[2] - μ1_r13*x_stm[2]       - μ2_r23*x_stm[2]
    dx_stm[6] =                        - μ1_r13*x_stm[3]       - μ2_r23*x_stm[3]
    # derivatives of STM
    A = dfdx_cr3bp(x_stm[1:6], μ, t)
    dx_stm[7:42] = reshape((A * reshape(x_stm[7:42],6,6)')', 36)
    return nothing
end


"""
    eom!(dx, x, p::CR3BP, t)

EOM for CR3BP system compatible with `DifferentialEquations.jl` in mutating form.
"""
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


"""
    dfdx(x, p::CR3BP, t)

Compute Jacobian of CR3BP EOM.
"""
function dfdx(x, p::CR3BP, t)
    A = zeros(6,6)
    A[1:3,4:6] = I(3)
    A[4,5] = 2.0
    A[5,4] = -2.0
    G1 = gravity_gradient_matrix(1-p.μ, [x[1] + p.μ, x[2], x[3]])
    G2 = gravity_gradient_matrix(p.μ, [x[1] - 1 + p.μ, x[2], x[3]])
    A[4:6,1:3] = [1 0 0; 0 1 0; 0 0 0] + G1 + G2
    return A
end


"""
    eom_stm!(dx_stm, x_stm, p::CR3BP, t)

EOM for CR3BP system with STM compatible with `DifferentialEquations.jl` in mutating form.
The state `x_stm` consists of the Cartesian r, v vectors & STM flattened in row-major order.
"""
function eom_stm!(dx_stm, x_stm, p::CR3BP, t)
    # calcualte intermediate scalar quantities
    r1 = sqrt((x_stm[1] + p.μ)^2 + x_stm[2]^2 + x_stm[3]^2)
    r2 = sqrt((x_stm[1] - 1 + p.μ)^2 + x_stm[2]^2 + x_stm[3]^2)
    μ1_r13 = (1 - p.μ) / r1^3
    μ2_r23 = p.μ / r2^3
    # derivatives of state
    dx_stm[1:3] = x_stm[4:6]
    dx_stm[4] =  2*x_stm[5] + x_stm[1] - μ1_r13*(x_stm[1] + p.μ) - μ2_r23*(x_stm[1] - 1 + p.μ)
    dx_stm[5] = -2*x_stm[4] + x_stm[2] - μ1_r13*x_stm[2]         - μ2_r23*x_stm[2]
    dx_stm[6] =                        - μ1_r13*x_stm[3]         - μ2_r23*x_stm[3]
    # derivatives of STM
    A = dfdx(x_stm[1:6], p, t)
    dx_stm[7:42] = reshape((A * reshape(x_stm[7:42],6,6)')', 36)
    return nothing
end