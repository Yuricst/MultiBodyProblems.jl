"""Utility functions"""

function outer(a::Vector, b::Vector)
    return a.*b'
end


function gravity_gradient_matrix(μ::Real, r::Vector)
    rnorm = norm(r)
    return μ/rnorm^5 * (3*outer(r, r) - I(3)*rnorm^2)
end