module PseudoSpectralReactionDiffusion
export PseudoSpectralProblem, PseudoSpectralSolution, PseudoSpectralIntegrator, init, steady_state_callback, x, step!, step_to!, get_sol, get_u

import Base: getindex, eachindex, lastindex
export getindex, eachindex, lastindex

import SciMLBase: EnsembleProblem, solve, remake, successful_retcode, DEIntegrator
export EnsembleProblem, solve, remake, successful_retcode

import SciMLBase
using SciMLBase: SplitODEProblem, ODEProblem, ODESolution, ODEFunction, update_coefficients!, ReturnCode, DiscreteCallback, terminate!, get_du
using SciMLBase.ReturnCode: Terminated
using SciMLOperators: DiagonalOperator
using OrdinaryDiffEqExponentialRK: ETDRK4
using FFTW: plan_r2r!, REDFT00, MEASURE, ScaledPlan
using Symbolics: variable, @variables, Num, sparsejacobian, build_function, substitute, get_variables
using Random: default_rng, AbstractRNG

"Spatial variable x∈[0,1]."
const x = variable(:x) |> Num

 # BC is Nothing for homogeneous for BC or Function for Heterogeneous BC.
 mutable struct PseudoSpectralProblem
    ode_problem::ODEProblem
    dims::Tuple{Int,Int}
    species::Vector{Num}
    reaction_params::Vector{Num}
    diffusion_params::Vector{Num}
    boundary_params::Vector{Num}
    initial_params::Vector{Num}
    plan::ScaledPlan
    initial_function::Function
    lifting_function::Union{Nothing, Function}
    rng::AbstractRNG
end

"""
    PseudoSpectralSolution
Solution object for PsuedoSpectralProblem.

# Indexing
- By time-step `sol[3]`.
- By species `sol[U]`.
"""
struct PseudoSpectralSolution
    species::Vector{Num}
    u::Vector{Matrix{Float64}}
    x:: Vector{Float64}
    t::Vector{Float64}
    retcode::ReturnCode.T
end

struct Parameters
    u :: Matrix{Float64} # Working array for dct.
    r :: Vector{Float64} # Reaction parameters.
    d :: Vector{Float64} # Diffusion parameters.
    ϕ :: Matrix{Float64} # Boundary lifting function
    Δϕ :: Matrix{Float64}
end


"""
    PseudoSpectralProblem(species, reaction_rates, diffusion_rates, boundary_conditions, initial_conditions, num_verts; p=nothing, noise=1e-4, rng=default_rng(), kwargs...)
Construct a PsuedoSpectralProblem object representing a reaction diffusion system of the form uₓₓ(x,t) = Duₜ(x,t) + f(u(x,t)).

# Arguments
PseudoSpectral expects Symbolics.jl expressions as inputs. The special variable 'x' ∈ [0,1] represents the spatial coordinate. Any variables other than 'x' and those supplied in `species` will be interpreted as parameters.

- `species`: Vector of variables corresponding to the components of u.
- `reaction_rates`: Vector of expressions representing f(u).
- `diffusion_rates`: Vector of Expressions representing diag(D).
- `boundary_conditions`: 2xn matrix of expressions representing Neumann boundary conditions. The two rows correspond to uₓ at the left and right boundaries.
- `initial_conditions`: Vector of expressions representing u(x,0).
- `num_verts`: Number of points in spatial discretisation.
- `p=nothing`: Dictionary associating parameters with numerical values.
- `noise=1e-4`: Guassian noise with σ²=`noise` is added to the initial conditions.
- `rng=default_rng()`: Random number generator for noise.
- `kwargs...`: Keyword arguments passed on to SciML's `solve`. For details see https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/.
"""
function PseudoSpectralProblem(species, reaction_rates, diffusion_rates, boundary_conditions, initial_conditions, num_verts; p=nothing, dt=0.1, noise=1e-4, rng=default_rng(), kwargs...)
    n = num_verts
    m = length(species)
    
    # Collect parameter symbols. 
    rs, ds, bs, is = (setdiff(collect_variables(exprs), x, species) for exprs in (reaction_rates, diffusion_rates, vec(boundary_conditions), initial_conditions))
    p = something(p, Dict(q => 0 for q in union(rs,ds,bs,is)))

    u = Matrix{Float64}(undef, n, m)
    plan = 1/sqrt(2*(n-1)) * plan_r2r!(u, REDFT00, 1; flags=MEASURE)

    fu0 = make_initial_function(initial_conditions, is, noise, n)
    lf = make_lifting_function(boundary_conditions, diffusion_rates, bs,ds, n)
    
    R = reaction_operator(species, reaction_rates, rs, plan, Val(!isnothing(lf)))
    D = diffusion_operator(diffusion_rates, ds, n)
    odeprob = SplitODEProblem(D, R, vec(u), Inf, nothing; dt, kwargs...)
    prob = PseudoSpectralProblem(odeprob, (n,m), species, rs, ds, bs, is, plan, fu0, lf, rng)
    remake(prob; p)
end

"""
    remake(prob::PseudoSpectralProblem; p=nothing, rng=nothing, kwargs...)

Return a new problem with updated parameters, random number generator, and/or solver options.
"""
function remake(prob::PseudoSpectralProblem; p=nothing, rng=nothing, kwargs...)
    if isnothing(p)
        prob.ode_problem = remake(prob.ode_problem; kwargs...)
        return prob
    end
    if !isnothing(rng)
        prob.rng=rng
    end
    r = Float64[p[k] for k in prob.reaction_params]
    d = Float64[p[k] for k in prob.diffusion_params]
    b = Float64[p[k] for k in prob.boundary_params]
    i = Float64[p[k] for k in prob.initial_params]

    w = Matrix{Float64}(undef,prob.dims...) # Allocate working memory for FFTW.
    u0 = prob.initial_function(i,prob.rng)
    lf = prob.lifting_function
    if !isnothing(lf)
        ϕ, Δϕ = lf(d,b)
        u0 .-= ϕ
    else
        ϕ = Δϕ = Matrix{Float64}(undef,0,0)
    end
    p = Parameters(w,r,d,ϕ,Δϕ)
    prob.plan * u0
    u0 = vec(u0)
    update_coefficients!(prob.ode_problem.f.f1.f, nothing, p, nothing) # Set parameter values in diffusion operator.
    prob.ode_problem = remake(prob.ode_problem; u0, p, kwargs...) # Set parameter values in SplitODEProblem.
    prob
end

"""
    solve(prob::PseudoSpectralProblem, alg=ETDRK4(); kwargs...)

See https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/.
Algorithm defaults to [ETDRK4](https://docs.sciml.ai/DiffEqDocs/stable/api/ordinarydiffeq/semilinear/ExponentialRK/#OrdinaryDiffEqExponentialRK.ETDRK4).
"""
function solve(prob::PseudoSpectralProblem, alg=ETDRK4(); kwargs...)
    odesol = SciMLBase.solve(prob.ode_problem, alg; kwargs...)
    PseudoSpectralSolution(prob, odesol)
end

# Separate constructor so we can use it both with solve and as an output function for EnsmbleProblem.
function PseudoSpectralSolution(prob::PseudoSpectralProblem, sol::ODESolution)
    u = [transform(prob, u) for u in sol.u]
    PseudoSpectralSolution(prob.species, u, range(0.0,1.0,prob.dims[1]), sol.t, sol.retcode)
end

function transform(prob::PseudoSpectralProblem, u)
    u = reshape(u, prob.dims)
    prob.plan * u
    if !isnothing(prob.lifting_function)
        u .+= prob.ode_problem.p.ϕ
    end
    u
end


function make_lifting_function(boundary_conditions, diffusion_rates, boundary_params,diffusion_params, n)
    iszero(boundary_conditions) && return nothing
    a,b = eachrow(boundary_conditions)
    X = range(0.0,1.0,n)
    ϕ = X.^2 * (b'-a')/2 + X * a'
    Δϕ = ((b-a).*diffusion_rates)'
    fϕ,_ = build_function(ϕ, boundary_params; expression=Val{false})
    fΔϕ,_ = build_function(Δϕ, diffusion_params, boundary_params; expression=Val{false})
    (d, b) -> (fϕ(b), fΔϕ(d,b))
end


function make_initial_function(initial_conditions, initial_params, initial_noise, n)
    m = length(initial_conditions)
    u0 = [substitute(ic, x=>X) for X in range(0,1,n), ic in initial_conditions]
    f,_= build_function(u0, initial_params; expression=Val{false})
    function (p,rng)
        noise = initial_noise * abs.(randn(rng,n,m))
        f(p) + noise
    end
end


"Build function for the reaction component, with `f(v+ϕ) + Δϕ` offset for non-zero-flux BCs."
function reaction_operator(species, reaction_rates, rs, plan!, ::Val{BC}) where BC
    n,m = size(plan!)
    @variables u[1:n, 1:m]
    # TODO: Clever things to make only spatially varying parameters expand?
    # Build an nxm matrix of derivatives, substituting reactants for u[i,j] and parameters for p[k,l].
    du = [substitute(expr, Dict([x=>X, zip(species,v)...])) for (v,X) in zip(eachrow(collect(u)), range(0,1,n)), expr in reaction_rates]
    _, f! = build_function(du, u, rs; expression=Val{false})
    
    function f̂!(du,u,p,t)
        du = reshape(du,n,m)
        copyto!(p.u, u)
        plan! * p.u
        BC && (p.u .+= p.ϕ)
        f!(du, p.u, p.r)
        BC && (du .+= p.Δϕ)
        plan! * du
        nothing
    end
    ODEFunction(f̂!)
end

"Build linear operator for the diffusion component."
function diffusion_operator(diffusion_rates, ps, n)
    k = 0:n-1 # Wavenumbers
    h = 1/(n-1)

    ## 2nd order Fourier differentiation coefficients.
    # For a continuous FT this would be σ² = (pi * k)^2, but corrected for
    # the discrete transform this becomes:
    σ² = @. ((2/h) * sin((h/2)*pi*k))^2

    λ = vec(-σ² * diffusion_rates') |> collect
    (f,f!) = build_function(λ, ps; expression=Val{false})
    λ0 = similar(λ, Float64)
    update!(λ,u,p,t) = f!(λ, p.d)
    DiagonalOperator(λ0; update_func! = update!)
end


function getindex(sol::PseudoSpectralSolution, species::Num)
    name=_nameof(species)
    i = findfirst(s -> _nameof(s.val)===name, sol.species)
    [vec(u[:,i]) for u in sol.u]
end

getindex(sol::PseudoSpectralSolution) = getindex(sol.u)
getindex(sol::PseudoSpectralSolution, i::Union{Int,CartesianIndex{1}}) = sol.u[i]

eachindex(sol::PseudoSpectralSolution) = eachindex(sol.u)
# (sol::PseudoSpectralSolution)(t) = transform(sol, sol.sol(t))

lastindex(sol::PseudoSpectralSolution) = lastindex(sol.u)

successful_retcode(sol::PseudoSpectralSolution) = SciMLBase.successful_retcode(sol.retcode)


"""
    EnsembleProblem(prob::PseudoSpectralProblem, params; output_func=nothing)

Construct an ensemble problem to solve the system in parallel for each of the supplied parameter sets.
"""
function EnsembleProblem(prob::PseudoSpectralProblem, params; output_func=nothing)
    prob_func(_prob,ctx) = remake(_prob; p=params[ctx.sim_id], rng=ctx.rng)
    EnsembleProblem(prob; prob_func, output_func)
end

"""
    EnsembleProblem(prob::PseudoSpectralProblem; prob_func, output_func=nothing)

Construct an ensemble problem to run the solver in parallel.
For details see https://docs.sciml.ai/DiffEqDocs/stable/features/ensemble/.
"""
function EnsembleProblem(prob::PseudoSpectralProblem; prob_func, output_func=nothing)
    _prob_func(_prob, ctx) = prob_func(prob, ctx).ode_problem
    function _output_func(sol, ctx) 
        ps_sol = PseudoSpectralSolution(prob, sol)
        isnothing(output_func) ?  (ps_sol,false) : output_func(ps_sol,ctx)
    end
    SciMLBase.EnsembleProblem(prob.ode_problem; prob_func=_prob_func, output_func=_output_func)
end


## Integrator interface
mutable struct PseudoSpectralIntegrator
    integrator::DEIntegrator
    prob::PseudoSpectralProblem
    ss::Float64
end

"""
    PseudoSpectralIntegrator(prob::PseudoSpectralProblem; alg=ETDRK4(), kwargs...)
Initialize an integrator for the problem.
See https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/.
Algorithm defaults to [ETDRK4](https://docs.sciml.ai/DiffEqDocs/stable/api/ordinarydiffeq/semilinear/ExponentialRK/#OrdinaryDiffEqExponentialRK.ETDRK4).
"""
function init(prob::PseudoSpectralProblem; alg=ETDRK4(), kwargs...)
    prob = remake(prob; kwargs...)
    integrator = SciMLBase.init(prob.ode_problem, alg)
    PseudoSpectralIntegrator(integrator, prob, Inf)
end

"""
    get_sol(integrator::PseudoSpectralIntegrator)
Return a solution object for the current integrator state.
"""
get_sol(integrator::PseudoSpectralIntegrator) = PseudoSpectralSolution(integrator.prob, integrator.integrator.sol)

"""
    get_u(integrator::PseudoSpectralIntegrator)
Return solution values at time `t`, stepping the integrator as necessary.
"""
function get_u(integrator::PseudoSpectralIntegrator, t) 
    step_to!(integrator, t)
    u = integrator.integrator.sol(t)
    transform(integrator.prob, u)
end

"""
    step!(integrator::PseudoSpectralIntegrator, dt=nothing, stop_at_tdt=false)
Advance the iterator by `dt`.
"""
function step!(integrator::PseudoSpectralIntegrator, dt=nothing, stop_at_tdt=false)
    SciMLBase.step!(integrator.integrator, dt, stop_at_tdt)
    if integrator.integrator.sol.retcode == Terminated
        integrator.ss = integrator.integrator.t
    end
end

"""
    step_to!(integrator::PseudoSpectralIntegrator, t, stop_at_tdt=false)
Advance the iterator to time `t`.
"""
function step_to!(integrator::PseudoSpectralIntegrator, t, stop_at_tdt=false)
    dt = max(0.0, t - integrator.integrator.t)
    step!(integrator, dt, stop_at_tdt)
end

"""
    remake(integrator::PseudoSpectralIntegrator; kwargs...)
Return a new integrator updated with `remake(integrator.prob; kwargs...)`.
"""
function remake(integrator::PseudoSpectralIntegrator; kwargs...)
    prob = remake(integrator.prob; kwargs...)
    init(prob; alg=integrator.integrator.alg)
end

"""
    steady_state_callback(tol=1e-4)
Callback function to be passed to `solve` to detect steady state. Terminates solver when |uₜ| ≤ `tol`.
"""
function steady_state_callback(tol=1e-4)
    condition(u,t,integrator) = isapprox(get_du(integrator), zero(u); atol=tol)
    DiscreteCallback(condition, terminate!)
end


## Symbolics utility functions
"Sort parameters by name."
sort_variables(p) = sort(p, by=_nameof)
#_nameof(v) = isspecies(v) ? nameof(v.f) : nameof(v)
function _nameof(v)  # TODO: Something less hacky.
    try
        nameof(v.f)
    catch e
        try
            nameof(v)
        catch
            nameof(v.val.f)
        end
    end
end


"Extract variables from a (possibly nested) collection of expressions and sort them by name."
collect_variables(exprs...) = collect_variables(exprs) # Combine multiple arguments.
collect_variables(exprs::Union{Tuple,Vector}) = exprs .|> collect_variables |> splat(union) |> sort_variables
collect_variables(expr) = get_variables(expr) |> collect # Call recursively until we get down to a single expression.

end

