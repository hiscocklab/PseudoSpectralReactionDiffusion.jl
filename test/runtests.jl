module PseudoSpectralTest

using PseudoSpectralReactionDiffusion
using Symbolics: @variables
using OrdinaryDiffEqExponentialRK: ETDRK4
using Test


@testset "zero flux" begin
    @variables U,g0,g1,d,a,b
    R = [0]
    D = [1/(pi)^2] # Divide by pi^2 for a domain of size pi.
    n=128
    dt=0.001
    B = [0,0]
    IC = [cos(pi*x)]
    prob = PseudoSpectralProblem([U], R, D, B, IC, n)
    sol = solve(prob, ETDRK4(); tspan=(0.0,2.0), dt=dt)
    @test successful_retcode(sol)
    @test sol[U][end] ≈ exp(-sol.t[end])*cos.(pi*sol.x) rtol=1e-2;
end
@testset "non-zero flux" begin
    @variables U,g0,g1,d,a,b
    R = [0]
    D = [1/(pi)^2] # Divide by pi^2 for a domain of size pi.
    n=128
    dt=0.001
    B = [pi,pi]
    IC = [pi*x]
    prob = PseudoSpectralProblem([U], R, D, B, IC, n)
    sol = solve(prob, ETDRK4(); tspan=(0.0,2.0), dt=dt)
    @test successful_retcode(sol)
    @test sol[U][end] ≈ (pi*sol.x) rtol=1e-2;
end
@testset "non-negative" begin
    @variables U,g0,g1,d,a,b
    R = [0]
    D = [1/(pi)^2] # Divide by pi^2 for a domain of size pi.
    n=128
    dt=0.001
    B = [-pi,0] # Inward flux
    IC = [1.0]
    prob = PseudoSpectralProblem([U], R, D, B, IC, n)
    sol = solve(prob, ETDRK4(); tspan=(0.0,2.0), dt=dt)
    @test successful_retcode(sol)
    @test all(>(0), sol[end])
end
@testset "EnsembleProblem" begin
    @variables U,g0,g1,d,a,b
    R = [0]
    D = [1/(pi)^2] # Divide by pi^2 for a domain of size pi.
    n=128
    dt=0.001
    @variables D
    B = [0,0]
    IC = [cos(pi*x)]
    prob = PseudoSpectralProblem([U], R, [D], B, IC, n; noise=0.0)
    # output_func(sol,ctx) = sol[U]
    params = [Dict(D=>d) for d in (1:3)/pi^2]  # Divide by pi^2 for a domain of size pi.
    ensembleprob = EnsembleProblem(prob, params)
    sol = solve(ensembleprob, ETDRK4(); tspan=(0.0,2.0), dt=dt, trajectories=length(params))
    sol1=sol.u[1]
    @test successful_retcode(sol1)
    @test sol1.u[end] ≈ exp(-sol1.t[end])*cos.(pi*sol1.x) rtol=1e-2;
end
@testset "Integrator" begin
    @variables U,g0,g1,d,a,b
    R = [0]
    D = [1/(pi)^2] # Divide by pi^2 for a domain of size pi.
    n=128
    dt=0.001
    R = [0]
    D = [1/(pi)^2] # Divide by pi^2 for a domain of size pi.
    n=128
    dt=0.001
    B = [0,0]
    IC = [cos(pi*x)]
    prob = PseudoSpectralProblem([U], R, D, B, IC, n; noise=0.0)
    int = init(prob)
    step_to!(int, 2.0)
    u = get_u(int,2.0)
    X=range(0.0,1.0,n)
    @test u ≈ exp(-2.0)*cos.(pi*X) rtol=1e-2;
end


end;