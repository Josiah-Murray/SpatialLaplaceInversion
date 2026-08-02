# SpatialLaplaceInversion.jl
NILaplace.jl provides a selection of algorithms for the inversion of the Laplace transform. It has a particular focus on applications with both temporal and spatial dimensions, but will perform equally well on time-only problems with no loss in efficiency.

The following algorithms are provided:
* **Weeks' method**: Expands the time-domain solution in terms of Laguerre polynomials. Efficient when the solution is needed at many times, but may fail to converge for challenging problems.
* **Gaver-Wynn rho (GWR)**: Approximates the time-domain solution at a particular time using a sequence derived by Gaver and based on Post's inversion formula. A sequence accelerator (the Wynn rho algorithm) is then applied. This is a powerful method and of the algorithms in this package, it generally gives the best convergence result. However, it can be computationally expensive and requires the use of arbitrary precision for large parameter values.
* **Direct quadrature method**: Approximates the time-domain solution by evaluating the Bromwich integral as a Riemann sum. This is a simple method with slow convergence relative to computational cost, however, it is useful for comparisons and may still converge when other methods (i.e. Weeks') fail. Suffers from Gibbs phenomenon when $f(0^+)\ne 0$.

## Installation

SpatialLaplaceInversion.jl can be installed using the command
```julia
pkg> add https://github.com/Josiah-Murray/SpatialLaplaceInversion.git
```

## Picking up from InverseLaplace.jl
I am indebted to the makers of the InverseLaplace.jl package for their work. The Gaver-Wynn rho algorithm provided here is a re-write based on their algorithm alongside the original papers by Abate and Valkó [1,2], and my implementation of Weeks' method also drew upon on theirs.

With my thanks conveyed, why am I making a new package when InverseLaplace.jl already exists? There are two main reasons:
1. The algorithms in use were not immeadiately convenient for use in problems with a spatial component (hence the new package's name).
2. Several of the features were only partially implemented or had small technical errors and (at time of writing) InverseLaplace.jl was not being maintained.

Whilst these were the motivating factors for development, I have added several features which set SpatialLaplaceInversion.jl apart from its predecessor:
* Capacity of all algorithms to handle vector valued Laplace domain functions i.e. $f(t)\in\mathbb R ^n$, $\bar f(s)\in\mathbb C^n$. This makes them convenient for use in space-time problems.
* The Gaver-Wynn rho algroithm allows for non-standard input types so that specialised multi-precision types can be used. (It is worth noting that the level of precision must be set *globally* for Julia's BigFloat data type. If working in multiple levels of precision, one would have to use, say, ArbNumerics.jl. This package allows for that.)
* The Gaver-Wynn rho algorithm has a small technical fault corrected (though this only marginally improves performance).
* A full suite of transform pairs are made available for testing purposes, based on the tables in 'Numerical methods for Laplace transform inversion' by Alan M. Cohen.

[1] Abate, J., & Valkó, P. P. (2004).
Multi‐precision Laplace transform inversion.
International Journal for Numerical Methods in Engineering, 60(5), 979–993.
https://doi.org/10.1002/nme.995

[2] Valkó, P. P., & Abate, J. (2004).
Comparison of sequence accelerators for the Gaver method of numerical Laplace transform inversion.
Computers & Mathematics with Applications, 48(3–4), 629–636.
https://doi.org/10.1016/j.camwa.2002.10.017
