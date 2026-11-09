using PseudoSpectralReactionDiffusion
using Symbolics: @variables
using OrdinaryDiffEqExponentialRK: ETDRK4
using Test

@testset "time-dependent reactions" begin
    @variables U a

    prob = PseudoSpectralProblem(
        [U],
        [a * t],
        [0],
        zeros(2, 1),
        [0],
        16;
        p = Dict(a => 2.0),
        noise = 0.0,
    )

    sol = solve(prob, ETDRK4(); tspan = (2.0, 3.0), dt = 0.001)

    @test successful_retcode(sol)
    @test sol[U][end] ≈ fill(5.0, 16) rtol = 1.0e-5
end
