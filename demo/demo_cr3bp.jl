"""Test for CR3BP system"""

using LinearAlgebra
using OrdinaryDiffEq
using GLMakie

using MultiBodyProblems

# define system
DU = 389703.0     # km
TU = 382981.0     # sec
μ = 1.215058560962404e-02
system = MultiBodyProblems.CR3BP(μ; DU=DU, TU=TU)

# propagate initial conditions
# JPL L4 Short orbit 1145
rv0 = [ 4.8784941344943100E-1,
		8.9915917622773556E-1,
		2.2507034007099641E-24,
		5.1985214035751473E-2,
		-2.1209779442896915E-2,
		-1.2508220242789951E-24]
period = 6.5821031915048422

# solve ODE over one period
prob = ODEProblem(MultiBodyProblems.eom!, rv0, (0.0, period), system)   # eom implied by parameters
sol = solve(prob, Vern7(), reltol=1e-13, abstol=1e-14)

# plot results
fig = Figure(size=(700,400))
ax = Axis3(fig[1,1], aspect=:data)
lines!(ax, Array(sol)[1,:], Array(sol)[2,:], Array(sol)[3,:])
display(fig)