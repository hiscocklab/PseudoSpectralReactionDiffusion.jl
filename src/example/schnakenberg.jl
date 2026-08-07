using PseudoSpectral
using Symbolics: @variables
using OrdinaryDiffEqExponentialRK: ETDRK4
using Plots

@variables U,V, Dᵤ, Dᵥ, a, b,L
R = [a-U+U^2*V, b-U^2*V]
D = [Dᵤ/L^2, Dᵥ/L^2]
B = [0 0; 0 0]
IC = [0,0]
n=128
dt=0.001

prob = PseudoSpectralProblem([U,V], R, D, B, IC, n; p = Dict(Dᵤ=>1.0, Dᵥ=>50.0, a=>0.2, b=>2.0, L=>50.0))
sol = solve(prob, ETDRK4(); tspan=(0.0,200.0), dt=dt)
plot(sol[U][end]; xlabel="x", ylabel="U", legend=false, title="Schnakenberg Pattern")
