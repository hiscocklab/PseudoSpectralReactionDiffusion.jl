
<a id='PseudoSpectral.jl'></a>

<a id='PseudoSpectral.jl-1'></a>

# PseudoSpectral.jl


<!– [![Build Status](https://github.com/twhiscock/ReactionDiffusion.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/twhiscock/ReactionDiffusion.jl/actions/workflows/CI.yml?query=branch%3Amaster) [![Latest Release (for users)](https://img.shields.io/badge/docs-stable-blue.svg)](https://hiscocklab.github.io/ReactionDiffusion.jl/stable) [![Master (for developers)](https://img.shields.io/badge/docs-dev-blue.svg)](https://hiscocklab.github.io/ReactionDiffusion.jl/dev) –>


This package provides a DCT based discretisation for 1D reaction-diffusion PDEs, designed to be used with the exponential time-differencing solvers provided by [OrdinaryDiffEqExponentialRK](https://docs.sciml.ai/DiffEqDocs/stable/api/ordinarydiffeq/semilinear/ExponentialRK/). It supports PDE systems of the form:


$ \begin{align*} \mathbf{u}*{t}(x,t) &= \mathbf{D} \mathbf{u}*{xx}(x,t) + \mathbf{f}(\mathbf{u}(x,t)) \
  \mathbf{u}(x,0) &= \mathbf{g}(x) \
  \mathbf{u}*x(0,t) &= a \
  \mathbf{u}*x(1,t) &= b \
\end{align*} $


For an interactive front-end which integrates with the Catalyst chemical modelling DSL, see [ReactionDiffusion.jl](https://github.com/hiscocklab/ReactionDiffusion.jl).


<a id='Example-Heat-Equation'></a>

<a id='Example-Heat-Equation-1'></a>

## Example - Heat Equation


PseudoSpectral.jl follows the interface of SciML:


```@example
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
```


<a id='API'></a>

<a id='API-1'></a>

## API

<a id='PseudoSpectral.PseudoSpectralProblem' href='#PseudoSpectral.PseudoSpectralProblem'>#</a>
**`PseudoSpectral.PseudoSpectralProblem`** &mdash; *Type*.



```julia
PseudoSpectralProblem(species, reaction_rates, diffusion_rates, boundary_conditions, initial_conditions, num_verts; p=nothing, noise=1e-4, rng=default_rng(), kwargs...)
```

Construct a PsuedoSpectralProblem object representing a reaction diffusion system of the form uₓₓ(x,t) = Duₜ(x,t) + f(u(x,t)).

**Arguments**

PseudoSpectral expects Symbolics.jl expressions as inputs. The special variable 'x' ∈ [0,1] represents the spatial coordinate. Any variables other than 'x' and those supplied in `species` will be interpreted as parameters.

  * `species`: Vector of variables corresponding to the components of u.
  * `reaction_rates`: Vector of expressions representing f(u).
  * `diffusion_rates`: Vector of Expressions representing diag(D).
  * `boundary_conditions`: 2xn matrix of expressions representing Neumann boundary conditions. The two rows correspond to uₓ at the left and right boundaries.
  * `initial_conditions`: Vector of expressions representing u(x,0).
  * `num_verts`: Number of points in spatial discretisation.
  * `p=nothing`: Dictionary associating parameters with numerical values.
  * `noise=1e-4`: Guassian noise with σ²=`noise` is added to the initial conditions.
  * `rng=default_rng()`: Random number generator for noise.
  * `kwargs...`: Keyword arguments passed on to SciML's `solve`. For details see https://docs.sciml.ai/DiffEqDocs/stable/basics/common*solver*opts/.


<a target='_blank' href='https://github.com/hiscocklab/PseudoSpectral/blob/bf6411a85421fcdd4cc0de1bcd3520635eaaf31c/src/PseudoSpectral.jl#L61-L78' class='documenter-source'>source</a><br>

<a id='PseudoSpectral.PseudoSpectralSolution' href='#PseudoSpectral.PseudoSpectralSolution'>#</a>
**`PseudoSpectral.PseudoSpectralSolution`** &mdash; *Type*.



```julia
PseudoSpectralSolution
```

Solution object for PsuedoSpectralProblem.

**Indexing**

  * By time-step `sol[3]`.
  * By species `sol[U]`.


<a target='_blank' href='https://github.com/hiscocklab/PseudoSpectral/blob/bf6411a85421fcdd4cc0de1bcd3520635eaaf31c/src/PseudoSpectral.jl#L36-L43' class='documenter-source'>source</a><br>

<a id='CommonSolve.solve' href='#CommonSolve.solve'>#</a>
**`CommonSolve.solve`** &mdash; *Function*.



```julia
solve(prob::PseudoSpectralProblem, alg=ETDRK4(); kwargs...)
```

See https://docs.sciml.ai/DiffEqDocs/stable/basics/common*solver*opts/. Algorithm defaults to [ETDRK4](https://docs.sciml.ai/DiffEqDocs/stable/api/ordinarydiffeq/semilinear/ExponentialRK/#OrdinaryDiffEqExponentialRK.ETDRK4).


<a target='_blank' href='https://github.com/hiscocklab/PseudoSpectral/blob/bf6411a85421fcdd4cc0de1bcd3520635eaaf31c/src/PseudoSpectral.jl#L135-L140' class='documenter-source'>source</a><br>

<a id='PseudoSpectral.PseudoSpectralIntegrator' href='#PseudoSpectral.PseudoSpectralIntegrator'>#</a>
**`PseudoSpectral.PseudoSpectralIntegrator`** &mdash; *Type*.



```julia
PseudoSpectralIntegrator(prob::PseudoSpectralProblem; alg=ETDRK4(), kwargs...)
```

Initialize an integrator for the problem. See https://docs.sciml.ai/DiffEqDocs/stable/basics/integrator/. Algorithm defaults to [ETDRK4](https://docs.sciml.ai/DiffEqDocs/stable/api/ordinarydiffeq/semilinear/ExponentialRK/#OrdinaryDiffEqExponentialRK.ETDRK4).


<a target='_blank' href='https://github.com/hiscocklab/PseudoSpectral/blob/bf6411a85421fcdd4cc0de1bcd3520635eaaf31c/src/PseudoSpectral.jl#L275-L280' class='documenter-source'>source</a><br>

<a id='PseudoSpectral.get_u' href='#PseudoSpectral.get_u'>#</a>
**`PseudoSpectral.get_u`** &mdash; *Function*.



```julia
get_u(integrator::PseudoSpectralIntegrator)
```

Return solution values at time `t`, stepping the integrator as necessary.


<a target='_blank' href='https://github.com/hiscocklab/PseudoSpectral/blob/bf6411a85421fcdd4cc0de1bcd3520635eaaf31c/src/PseudoSpectral.jl#L292-L295' class='documenter-source'>source</a><br>

<a id='PseudoSpectral.get_sol' href='#PseudoSpectral.get_sol'>#</a>
**`PseudoSpectral.get_sol`** &mdash; *Function*.



```julia
get_sol(integrator::PseudoSpectralIntegrator)
```

Return a solution object for the current integrator state.


<a target='_blank' href='https://github.com/hiscocklab/PseudoSpectral/blob/bf6411a85421fcdd4cc0de1bcd3520635eaaf31c/src/PseudoSpectral.jl#L286-L289' class='documenter-source'>source</a><br>

<a id='PseudoSpectral.step!' href='#PseudoSpectral.step!'>#</a>
**`PseudoSpectral.step!`** &mdash; *Function*.



```julia
step!(integrator::PseudoSpectralIntegrator, dt=nothing, stop_at_tdt=false)
```

Advance the iterator by `dt`.


<a target='_blank' href='https://github.com/hiscocklab/PseudoSpectral/blob/bf6411a85421fcdd4cc0de1bcd3520635eaaf31c/src/PseudoSpectral.jl#L302-L305' class='documenter-source'>source</a><br>

<a id='PseudoSpectral.step_to!' href='#PseudoSpectral.step_to!'>#</a>
**`PseudoSpectral.step_to!`** &mdash; *Function*.



```julia
step_to!(integrator::PseudoSpectralIntegrator, t, stop_at_tdt=false)
```

Advance the iterator to time `t`.


<a target='_blank' href='https://github.com/hiscocklab/PseudoSpectral/blob/bf6411a85421fcdd4cc0de1bcd3520635eaaf31c/src/PseudoSpectral.jl#L313-L316' class='documenter-source'>source</a><br>


<a id='Support,-citation-and-future-developments'></a>

<a id='Support,-citation-and-future-developments-1'></a>

## Support, citation and future developments


If you find ReactionDiffusion.jl helpful in your research, teaching, or other activities, please star the repository and consider citing [this paper](https://www.biorxiv.org/content/10.1101/2025.05.27.656324v1).


We are a small team of academic researchers from the [Hiscock Lab](https://twhiscock.github.io/), who build mathematical models of developing embryos and tissues. We have found these tools helpful in our own research, and make them available in case you find them helpful in your research too. We hope to extend the functionality of ReactionDiffusion.jl as our future projects, funding and time allows.


This work is supported by ERC grant SELFORG-101161207, and UK Research and Innovation (Biotechnology and Biological Sciences Research Council, grant number BB/W003619/1) 


*Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them*


![ERC_logo](assets/LOGO_ERC-FLAG_FP.png)

