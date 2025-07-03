module MultiBodyProblems

using DataFrames
using JSON
using LinearAlgebra
using OrdinaryDiffEq
using Printf

include("utils.jl")
include("CR3BP.jl")
include("BCR4BP.jl")
include("periodic_orbit.jl")


export CR3BP
export eom_cr3bp!, dfdx_cr3bp, eom_stm_cr3bp!

export BCR4BP
export eom_bcr4bp!, dfdx_bcr4bp, eom_stm_bcr4bp!

export eom!, eom_stm!   # generic eom function names accepting systems (e.g. CR3BP, BCR4BP) as parameters


end # module MultiBodyProblems
