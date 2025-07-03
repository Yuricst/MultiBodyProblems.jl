"""Bi-Circular Restricted Four-Body Problem"""

"""
    BCR4BP(μ1::Float64, μ2::Float64; DU::Float64 = 1.0, TU::Float64 = 1.0)

Struct for BCR4BP system.
"""
mutable struct BCR4BP <: MultiBodyProblem
    μ::Float64
    μ_s::Float64
    a_s::Float64
    ω_s::Float64

    μ1::Float64
    μ2::Float64
    radius_m1::Union{Nothing,Float64}
    radius_m2::Union{Nothing,Float64}
    DU::Float64
    TU::Float64
    VU::Float64

    function BCR4BP(
        μ::Float64,
        μ_s::Float64,
        a_s::Float64,
        ω_s::Float64;
        DU::Float64 = 1.0,
        TU::Float64 = 1.0,
        radius_m1::Union{Nothing,Float64} = nothing,
        radius_m2::Union{Nothing,Float64} = nothing
    )
        μ1 = 1 - μ
        μ2 = μ
        VU = DU / TU
        return new(μ, μ_s, a_s, ω_s, μ1, μ2, radius_m1, radius_m2, DU, TU, VU)
    end
end


function Base.show(io::IO, bcr4bp::BCR4BP)
    println(io, "BCR4BP system")
    @printf(io, "   μ         = %f\n", bcr4bp.μ)
    @printf(io, "   μ_s       = %f\n", bcr4bp.μ_s)
    @printf(io, "   a_s       = %f\n", bcr4bp.a_s)
    @printf(io, "   ω_s       = %f\n", bcr4bp.ω_s)
    if bcr4bp.radius_m1 !== nothing
        @printf(io, "   radius_m1 = %f\n", bcr4bp.radius_m1)
    end
    if bcr4bp.radius_m2 !== nothing
        @printf(io, "   radius_m2 = %f\n", bcr4bp.radius_m2)
    end
    @printf(io, "   DU        = %f\n", bcr4bp.DU)
    @printf(io, "   TU        = %f\n", bcr4bp.TU)
end


"""
    eom_bcr4bp!(dx, x, p::Vector, t)

EOM for BCR4BP system compatible with `DifferentialEquations.jl` in mutating form.
Parameters `p` is a vector containing: `[μ, μ_s, a_s, ω_s]`.
"""
function eom_bcr4bp!(dx, x, p::Vector, t)
    # unpack parameters
    μ, μ_s, a_s, ω_s = p

    # calculate radii to first and second primaries
    r1 = sqrt( (x[1]+μ)^2 + x[2]^2 + x[3]^2 )
    r2 = sqrt( (x[1]-1+μ)^2 + x[2]^2 + x[3]^2 )

    # sun position
    xs = a_s*cos(ω_s*t)
    ys = a_s*sin(ω_s*t)
    zs = 0.0
    rs = sqrt( (x[1]-xs)^2 + (x[2]-ys)^2 + (x[3]-zs)^2 )

    # intermediate scalar quantities
    μ1_r13 = (1 - μ) / r1^3
    μ2_r23 = μ / r2^3
    μs_rs3 = μ_s / rs^3
    μs_a3  = μ_s / a_s^3

    # state derivatives
    dx[1:3] = x[4:6]
    dx[4] =  2*x[5] + x[1] - μ1_r13*(μ+x[1]) - μ2_r23*(x[1]-1+μ) + ( -μs_rs3*(x[1]-xs) - μs_a3*xs )
    dx[5] = -2*x[4] + x[2] - μ1_r13*x[2]     - μ2_r23*x[2]       + ( -μs_rs3*(x[2]-ys) - μs_a3*ys )
    dx[6] =                - μ1_r13*x[3]     - μ2_r23*x[3]       + ( -μs_rs3*(x[3])    - μs_a3*zs )
    return nothing
end


"""
    dfdx_bcr4bp(x, p::Vector, t)

Compute Jacobian of BCR4BP EOM.
Parameters `p` is a vector containing: `[μ, μ_s, a_s, ω_s]`.
"""
function dfdx_bcr4bp(x, p::Vector, t)
    # unpack parameters
    μ, μ_s, a_s, ω_s = p
    A = zeros(6,6)
    A[1:3,4:6] = I(3)
    A[4,5] = 2.0
    A[5,4] = -2.0
    G1 = gravity_gradient_matrix(1-μ, [x[1] + μ, x[2], x[3]])
    G2 = gravity_gradient_matrix(μ, [x[1] - 1 + μ, x[2], x[3]])

    # sun position
    xs = a_s*cos(ω_s*t)
    ys = a_s*sin(ω_s*t)
    zs = 0.0
    rs = sqrt( (x[1]-xs)^2 + (x[2]-ys)^2 + (x[3]-zs)^2 )

    # intermediate scalar quantities
    μs_rs5 = μ_s / rs^5
    μs_rs3 = μ_s / rs^3

    # potential terms due to fourth body
    Ups_xx = 3*μs_rs5*(x[1]-xs)^2 - μs_rs3
    Ups_yy = 3*μs_rs5*(x[2]-ys)^2 - μs_rs3
    Ups_zz = 3*μs_rs5*(x[3]-zs)^2 - μs_rs3
    Ups_xy = 3*μs_rs5*(x[2]-ys)*(x[1]-xs)
    Ups_xz = 3*μs_rs5*(x[3]-zs)*(x[1]-xs)
    Ups_yz = 3*μs_rs5*(x[3]-zs)*(x[2]-ys)
    U_sun = [Ups_xx Ups_xy Ups_xz;
             Ups_xy Ups_yy Ups_yz;
             Ups_xz Ups_yz Ups_zz]
    A[4:6,1:3] = [1 0 0; 0 1 0; 0 0 0] + G1 + G2 + U_sun
    return A
end


"""
    eom_stm_bcr4bp!(dx_stm, x_stm, p, t)

EOM for BCR4BP system with STM compatible with `DifferentialEquations.jl` in mutating form.
Parameters `p` is a vector containing: `[μ, μ_s, a_s, ω_s]`.
The state `x_stm` consists of the Cartesian r, v vectors & STM flattened in row-major order.
"""
function eom_stm_bcr4bp!(dx_stm, x_stm, p::Vector, t)
    # unpack parameters
    μ, μ_s, a_s, ω_s = p

    # calculate radii to first and second primaries
    r1 = sqrt( (x_stm[1]+μ)^2 + x_stm[2]^2 + x_stm[3]^2 )
    r2 = sqrt( (x_stm[1]-1+μ)^2 + x_stm[2]^2 + x_stm[3]^2 )

    # sun position
    xs = a_s*cos(ω_s*t)
    ys = a_s*sin(ω_s*t)
    zs = 0.0
    rs = sqrt( (x_stm[1]-xs)^2 + (x_stm[2]-ys)^2 + (x_stm[3]-zs)^2 )

    # intermediate scalar quantities
    μ1_r13 = (1 - μ) / r1^3
    μ2_r23 = μ / r2^3
    μs_rs3 = μ_s / rs^3
    μs_a3  = μ_s / a_s^3
    μs_rs5 = μ_s / rs^5     # used in Jacobian

    # state derivatives
    dx_stm[1:3] = x_stm[4:6]
    dx_stm[4] =  2*x_stm[5] + x_stm[1] - μ1_r13*(μ+x_stm[1]) - μ2_r23*(x_stm[1]-1+μ) + ( -μs_rs3*(x_stm[1]-xs) - μs_a3*xs )
    dx_stm[5] = -2*x_stm[4] + x_stm[2] - μ1_r13*x_stm[2]     - μ2_r23*x_stm[2]       + ( -μs_rs3*(x_stm[2]-ys) - μs_a3*ys )
    dx_stm[6] =                        - μ1_r13*x_stm[3]     - μ2_r23*x_stm[3]       + ( -μs_rs3*(x_stm[3])    - μs_a3*zs )

    # derivatives of STM (we choose not to call dfdx_bcr4bp so we can re-use precomputed quantities)
    A = zeros(6,6)
    A[1:3,4:6] = I(3)
    A[4,5] = 2.0
    A[5,4] = -2.0
    G1 = gravity_gradient_matrix(1-μ, [x_stm[1] + μ, x_stm[2], x_stm[3]])
    G2 = gravity_gradient_matrix(μ, [x_stm[1] - 1 + μ, x_stm[2], x_stm[3]])

    # potential terms due to fourth body
    Ups_xx = 3*μs_rs5*(x_stm[1]-xs)^2 - μs_rs3
    Ups_yy = 3*μs_rs5*(x_stm[2]-ys)^2 - μs_rs3
    Ups_zz = 3*μs_rs5*(x_stm[3]-zs)^2 - μs_rs3
    Ups_xy = 3*μs_rs5*(x_stm[2]-ys)*(x_stm[1]-xs)
    Ups_xz = 3*μs_rs5*(x_stm[3]-zs)*(x_stm[1]-xs)
    Ups_yz = 3*μs_rs5*(x_stm[3]-zs)*(x_stm[2]-ys)
    U_sun = [Ups_xx Ups_xy Ups_xz;
             Ups_xy Ups_yy Ups_yz;
             Ups_xz Ups_yz Ups_zz]
    A[4:6,1:3] = [1 0 0; 0 1 0; 0 0 0] + G1 + G2 + U_sun

    dx_stm[7:42] = reshape((A * reshape(x_stm[7:42],6,6)')', 36)
    return nothing
end


"""
    eom!(dx, x, p::BCR4BP, t)

EOM for BCR4BP system compatible with `DifferentialEquations.jl` in mutating form.
"""
function eom!(dx, x, p::BCR4BP, t)
    # calculate radii to first and second primaries
    r1 = sqrt( (x[1]+p.μ)^2 + x[2]^2 + x[3]^2 )
    r2 = sqrt( (x[1]-1+p.μ)^2 + x[2]^2 + x[3]^2 )

    # sun position
    xs = p.a_s*cos(p.ω_s*t)
    ys = p.a_s*sin(p.ω_s*t)
    zs = 0.0
    rs = sqrt( (x[1]-xs)^2 + (x[2]-ys)^2 + (x[3]-zs)^2 )

    # intermediate scalar quantities
    μ1_r13 = (1 - p.μ) / r1^3
    μ2_r23 = p.μ / r2^3
    μs_rs3 = p.μ_s / rs^3
    μs_a3  = p.μ_s / p.a_s^3

    # state derivatives
    dx[1:3] = x[4:6]
    dx[4] =  2*x[5] + x[1] - μ1_r13*(p.μ+x[1]) - μ2_r23*(x[1]-1+p.μ) + ( -μs_rs3*(x[1]-xs) - μs_a3*xs )
    dx[5] = -2*x[4] + x[2] - μ1_r13*x[2]       - μ2_r23*x[2]         + ( -μs_rs3*(x[2]-ys) - μs_a3*ys )
    dx[6] =                - μ1_r13*x[3]       - μ2_r23*x[3]         + ( -μs_rs3*(x[3])    - μs_a3*zs )
    return nothing
end


"""
    dfdx(x, p::BCR4BP, t)

Compute Jacobian of BCR4BP EOM.
"""
function dfdx(x, p::BCR4BP, t)
    A = zeros(6,6)
    A[1:3,4:6] = I(3)
    A[4,5] =  2.0
    A[5,4] = -2.0
    G1 = gravity_gradient_matrix(1-p.μ, [x[1] + p.μ, x[2], x[3]])
    G2 = gravity_gradient_matrix(p.μ, [x[1] - 1 + p.μ, x[2], x[3]])

    # sun position
    xs = p.a_s*cos(p.ω_s*t)
    ys = p.a_s*sin(p.ω_s*t)
    zs = 0.0
    rs = sqrt( (x[1]-xs)^2 + (x[2]-ys)^2 + (x[3]-zs)^2 )

    # intermediate scalar quantities
    μs_rs5 = p.μ_s / rs^5
    μs_rs3 = p.μ_s / rs^3

    # potential terms due to fourth body
    Ups_xx = 3*μs_rs5*(x[1]-xs)^2 - μs_rs3
    Ups_yy = 3*μs_rs5*(x[2]-ys)^2 - μs_rs3
    Ups_zz = 3*μs_rs5*(x[3]-zs)^2 - μs_rs3
    Ups_xy = 3*μs_rs5*(x[2]-ys)*(x[1]-xs)
    Ups_xz = 3*μs_rs5*(x[3]-zs)*(x[1]-xs)
    Ups_yz = 3*μs_rs5*(x[3]-zs)*(x[2]-ys)
    U_sun = [Ups_xx Ups_xy Ups_xz;
             Ups_xy Ups_yy Ups_yz;
             Ups_xz Ups_yz Ups_zz]
    A[4:6,1:3] = [1 0 0; 0 1 0; 0 0 0] + G1 + G2 + U_sun
    return A
end


function eom_stm!(dx_stm, x_stm, p::BCR4BP, t)
    # calculate radii to first and second primaries
    r1 = sqrt( (x_stm[1]+p.μ)^2 + x_stm[2]^2 + x_stm[3]^2 )
    r2 = sqrt( (x_stm[1]-1+p.μ)^2 + x_stm[2]^2 + x_stm[3]^2 )

    # sun position
    xs = p.a_s*cos(p.ω_s*t)
    ys = p.a_s*sin(p.ω_s*t)
    zs = 0.0
    rs = sqrt( (x_stm[1]-xs)^2 + (x_stm[2]-ys)^2 + (x_stm[3]-zs)^2 )

    # intermediate scalar quantities
    μ1_r13 = (1 - p.μ) / r1^3
    μ2_r23 = p.μ / r2^3
    μs_rs3 = p.μ_s / rs^3
    μs_a3  = p.μ_s / p.a_s^3
    μs_rs5 = p.μ_s / rs^5     # used in Jacobian

    # state derivatives
    dx_stm[1:3] = x_stm[4:6]
    dx_stm[4] =  2*x_stm[5] + x_stm[1] - μ1_r13*(p.μ+x_stm[1]) - μ2_r23*(x_stm[1]-1+p.μ) + ( -μs_rs3*(x_stm[1]-xs) - μs_a3*xs )
    dx_stm[5] = -2*x_stm[4] + x_stm[2] - μ1_r13*x_stm[2]       - μ2_r23*x_stm[2]         + ( -μs_rs3*(x_stm[2]-ys) - μs_a3*ys )
    dx_stm[6] =                        - μ1_r13*x_stm[3]       - μ2_r23*x_stm[3]         + ( -μs_rs3*(x_stm[3])    - μs_a3*zs )

    # derivatives of STM (we choose not to call dfdx_bcr4bp so we can re-use precomputed quantities)
    A = zeros(6,6)
    A[1:3,4:6] = I(3)
    A[4,5] = 2.0
    A[5,4] = -2.0
    G1 = gravity_gradient_matrix(1-p.μ, [x_stm[1] + p.μ, x_stm[2], x_stm[3]])
    G2 = gravity_gradient_matrix(p.μ, [x_stm[1] - 1 + p.μ, x_stm[2], x_stm[3]])

    # potential terms due to fourth body
    Ups_xx = 3*μs_rs5*(x_stm[1]-xs)^2 - μs_rs3
    Ups_yy = 3*μs_rs5*(x_stm[2]-ys)^2 - μs_rs3
    Ups_zz = 3*μs_rs5*(x_stm[3]-zs)^2 - μs_rs3
    Ups_xy = 3*μs_rs5*(x_stm[2]-ys)*(x_stm[1]-xs)
    Ups_xz = 3*μs_rs5*(x_stm[3]-zs)*(x_stm[1]-xs)
    Ups_yz = 3*μs_rs5*(x_stm[3]-zs)*(x_stm[2]-ys)
    U_sun = [Ups_xx Ups_xy Ups_xz;
             Ups_xy Ups_yy Ups_yz;
             Ups_xz Ups_yz Ups_zz]
    A[4:6,1:3] = [1 0 0; 0 1 0; 0 0 0] + G1 + G2 + U_sun

    dx_stm[7:42] = reshape((A * reshape(x_stm[7:42],6,6)')', 36)
    return nothing
end


