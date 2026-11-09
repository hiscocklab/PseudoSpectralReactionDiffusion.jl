using PseudoSpectralReactionDiffusion
using Symbolics: @variables
using OrdinaryDiffEqExponentialRK: ETDRK4
using Test

@testset "zero flux" begin
    @variables U, g0, g1, d, a, b
    R = [0]
    D = [1 / (pi)^2] # Divide by pi^2 for a domain of size pi.
    n = 128
    dt = 0.001
    B = [0, 0]
    IC = [cos(pi * x)]
    prob = PseudoSpectralProblem([U], R, D, B, IC, n)
    sol = solve(prob, ETDRK4(); tspan = (0.0, 2.0), dt = dt)
    @test successful_retcode(sol)
    @test sol[U][end] ≈ exp(-sol.t[end]) * cos.(pi * sol.x) rtol = 1.0e-2
end

@testset "non-zero flux" begin
    @variables U, g0, g1, d, a, b
    R = [0]
    D = [1 / (pi)^2] # Divide by pi^2 for a domain of size pi.
    n = 128
    dt = 0.001
    B = [pi, pi]
    IC = [pi * x]
    prob = PseudoSpectralProblem([U], R, D, B, IC, n)
    sol = solve(prob, ETDRK4(); tspan = (0.0, 2.0), dt = dt)
    @test successful_retcode(sol)
    @test sol[U][end] ≈ (pi * sol.x) rtol = 1.0e-2
end

@testset "non-negative" begin
    @variables U, g0, g1, d, a, b
    R = [0]
    D = [1 / (pi)^2] # Divide by pi^2 for a domain of size pi.
    n = 128
    dt = 0.001
    B = [-pi, 0] # Inward flux
    IC = [1.0]
    prob = PseudoSpectralProblem([U], R, D, B, IC, n)
    sol = solve(prob, ETDRK4(); tspan = (0.0, 2.0), dt = dt)
    @test successful_retcode(sol)
    @test all(>(0), sol[end])
end

@testset "dealias" begin
    @variables U, g0, g1, d, a, b
    R = [0]
    D = [1 / (pi)^2] # Divide by pi^2 for a domain of size pi.
    n = 128
    dt = 0.001
    B = [0, 0]
    IC = [cos(pi * x)]
    prob = PseudoSpectralProblem([U], R, D, B, IC, n; dealias = true)
    sol = solve(prob, ETDRK4(); tspan = (0.0, 2.0), dt = dt)
    @test successful_retcode(sol)
    @test sol[U][end] ≈ exp(-sol.t[end]) * cos.(pi * sol.x) rtol = 1.0e-2
end
