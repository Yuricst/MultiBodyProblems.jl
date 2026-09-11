"""Ensemble CR3BP propagation of a bundle of initial states"""

using LinearAlgebra
using OrdinaryDiffEq
using Random
using GLMakie

using MultiBodyProblems

# define system
DU = 389703.0     # km
TU = 382981.0     # sec
μ = 1.215058560962404e-02
system = MultiBodyProblems.CR3BP(μ; DU=DU, TU=TU)

# nominal initial conditions
# JPL L4 Short orbit 1145
rv0 = [ 4.8784941344943100E-1,
		8.9915917622773556E-1,
		2.2507034007099641E-24,
		5.1985214035751473E-2,
		-2.1209779442896915E-2,
		-1.2508220242789951E-24]
period = 6.5821031915048422

# bundle of 100 initial states distributed around rv0
n_traj = 100
σ = 1e-4   # isotropic perturbation amplitude (nondimensional)
rng = Xoshiro(42)
u0s = [rv0 + σ * randn(rng, 6) for _ in 1:n_traj]

# template problem (same EOM, timespan, and system as the single-state demo)
prob = ODEProblem(MultiBodyProblems.eom!, rv0, (0.0, period), system)

function prob_func(prob, ctx)
    remake(prob, u0 = u0s[ctx.sim_id])
end

ensemble_prob = EnsembleProblem(prob; prob_func = prob_func)
sim = solve(ensemble_prob, Vern7(); trajectories = n_traj, reltol=1e-13, abstol=1e-14)

# plot the bundle
fig = Figure(size=(700,400))
ax = Axis3(fig[1,1], aspect=:data)
for soli in sim.u
    traj = Array(soli)
    lines!(ax, traj[1,:], traj[2,:], traj[3,:]; color=(:royalblue, 0.25))
end
display(fig)
