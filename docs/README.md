# PseudoSpectral.jl

<!-- [![Build Status](https://github.com/twhiscock/ReactionDiffusion.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/twhiscock/ReactionDiffusion.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Latest Release (for users)](https://img.shields.io/badge/docs-stable-blue.svg)](https://hiscocklab.github.io/ReactionDiffusion.jl/stable)
[![Master (for developers)](https://img.shields.io/badge/docs-dev-blue.svg)](https://hiscocklab.github.io/ReactionDiffusion.jl/dev) -->

This package provides a DCT based discretisation for 1D reaction-diffusion PDEs, designed to be used with the exponential time-differencing solvers provided by [OrdinaryDiffEqExponentialRK](https://docs.sciml.ai/DiffEqDocs/stable/api/ordinarydiffeq/semilinear/ExponentialRK/). It supports PDE systems of the form:

$ \begin{align*} \mathbf{u}_{t}(x,t) &= \mathbf{D} \mathbf{u}_{xx}(x,t) + \mathbf{f}(\mathbf{u}(x,t)) \\
  \mathbf{u}(x,0) &= \mathbf{g}(x) \\
  \mathbf{u}_x(0,t) &= a \\
  \mathbf{u}_x(1,t) &= b \\
\end{align*} $

For an interactive front-end which integrates with the Catalyst chemical modelling DSL, see [ReactionDiffusion.jl](https://github.com/hiscocklab/ReactionDiffusion.jl).

## Example - Heat Equation
PseudoSpectral.jl follows the interface of SciML:
```
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
![Schnakenberg system](schnakenberg.png)

## API

```@docs
PseudoSpectralProblem
PseudoSpectralSolution
solve
PseudoSpectralIntegrator
get_u
get_sol
step!
step_to!
```

## Support, citation and future developments

If you find ReactionDiffusion.jl helpful in your research, teaching, or other activities, please star the repository and consider citing [this paper](https://www.biorxiv.org/content/10.1101/2025.05.27.656324v1).

We are a small team of academic researchers from the [Hiscock Lab](https://twhiscock.github.io/), who build mathematical models of developing embryos and tissues. We have found these tools helpful in our own research, and make them available in case you find them helpful in your research too. We hope to extend the functionality of ReactionDiffusion.jl as our future projects, funding and time allows.

This work is supported by ERC grant SELFORG-101161207, and UK Research and Innovation (Biotechnology and Biological Sciences Research Council, grant number BB/W003619/1) 

*Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them*

![ERC_logo](docs/src/assets/LOGO_ERC-FLAG_FP.png)

