module MultiBodyProblems

using OrdinaryDiffEq
using Printf

include("CR3BP.jl")

export CR3BP
export eom_cr3bp!, dfdx_cr3bp, eom_stm_cr3bp!
export eom!, eom_stm!


end # module MultiBodyProblems
