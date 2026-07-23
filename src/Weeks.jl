"""
An implementation of Weeks method for numerical Laplace inversion. Weeks method approximates
the inversion of the Laplace domain function using a series of Laguerre
polynomials. Unlike the other methods, this results in an approximating
function that can be re-sampled at the required time values. The idea here is
that the coefficients of the approximating series are calculated once and stored in
a struct, which can then be passed into a separate function that evaluates the
approximation for a particular value of t.

The work flow for the method is as follows:
1. Create an instance of the mutable struct [`WeeksApproximation`](@ref) using the
function [`GenerateWeeksApproximation`](@ref) which computes the coefficients of the
approximating Laguerre series.

2. Approximate the inversion at a time t by calling `EvalWeeks` and
passing in the created instance of the [`WeeksApproximation`](@ref) struct, and
a time value `t`.

More details on the particular method used here can be found in:

Lyness, J. N., & Giunta, G. (1986).
A Modification of the Weeks Method for
Numerical Inversion of the Laplace Transform.
Mathematics of Computation, 47(175), 313.
https://doi.org/10.2307/2008097


Clenshaw, C. W. (1955). A
note on the summation of Chebyshev series.
Mathematics of Computation, 9(51), 118–120.
https://doi.org/10.1090/S0025-5718-1955-0071856-0

"""
module Weeks

"""
    WeeksApproximation(lapFunc, N, σ, b, coefficients, evalType)

Store the coefficients of a Weeks method approximation of order `N` to the Laplace domain function `lapFunc`, using parameters `σ` and `b`.
Computations are performed by converting numbers to type `evalType` (in general this will be `Float64` or `BigFloat`).
Should be initialised using the related function [`GenerateWeeksApproximation`](@ref).
"""
mutable struct WeeksApproximation
  lapFunc #Laplace domain function to be approximated
  N #Order of the leading term in the polynomial approximation.
  σ #Control parameter
  b #Control parameter
  coefficients #Coefficients of the Laguerre series approximating lapFunc
  evalType #The type in which to do the calculations. In general, use Float64. Also built to use BigFloat.
end

#MARK: Main functions

"""
    GenerateWeeksApproximation(lapFunc, N, σ, b; evalType = Float64)

Calculate the coefficients in the `N`th order Laguerre series of Weeks method for a Laplace domain function `lapFunc(s)`, using parameters `σ` and `b`,  and return a [`WeeksApproximation`](@ref) struct. By default, performs calculations in `Float64`, but other data types (e.g. `BigFloat`) can be used by changing the `evalType` keyword argument.

The inversion can be performed for a chosen time using [`EvalWeeks`](@ref).

# Examples

```jldoctest
julia> f = s -> 1/s^2
Julia> WeeksApprox = GenerateWeeksApproximation(f, 40, 0.5, 1.0)
julia> EvalWeeks(WeeksApprox, 2.0) #Evaluate at t=2.0
2.0000000000000018
```


```jldoctest
julia> f = s -> [1/s^2, 1/s^3]
Julia> WeeksApprox = GenerateWeeksApproximation(f, 40, 0.5, 1.0)
julia> EvalWeeks(WeeksApprox, 2.0) #Evaluate at t=2.0
2-element Vector{Float64}:
 2.0000000000000018
 2.000000000000003
```
"""
function GenerateWeeksApproximation(lapFunc, N, σ, b; evalType = Float64)
  σ = convert(evalType, σ)
  b = convert(evalType, b)
  coefficients = ComputeCoefficients(lapFunc, N, σ, b)
  return WeeksApproximation(lapFunc, N, σ, b, coefficients, evalType)
end



"""
    EvalWeeks(Weeks::WeeksApproximation, t)

Evaluate the approximation through Weeks method stored in `Weeks` (created using [`GenerateWeeksApproximation`](@ref)) at the time `t`, using the Clenshaw algorithm.

# Examples

```jldoctest
julia> f = s -> 1/s^2
Julia> WeeksApprox = GenerateWeeksApproximation(f, 40, 0.5, 1.0)
julia> EvalWeeks(WeeksApprox, 2.0) #Evaluate at t=2.0
2.0000000000000018
```
"""
function EvalWeeks(Weeks::WeeksApproximation, t)
  evalType = Weeks.evalType
  N = Weeks.N
  α_n = n -> - (1/(convert(evalType,n)+1))*t + (2*convert(evalType,n) + 1)/(convert(evalType, n)+1)
  β_n = n -> - convert(evalType, n)/(convert(evalType, n)+1)
  p_0 = 1
  p_1 = 1-t

  #Set initial terms in the Clenshaw algorithm.
  bp2 = 0
  bp1 = 0
  b = 0

  #Back propagation for Clenshaw algorithm
  for i in N:-1:1
    b = Weeks.coefficients[i+1] .+ α_n(i)*bp1 .+ β_n(i+1)*bp2
    bp2 = bp1
    bp1 = b

  end

  approximation = Weeks.coefficients[1] .+ (1-t)* bp1 .+ β_n(1)*bp2



  return approximation
end



#MARK: Internal functions

"""
    ComputeCoefficients(lapFunc, N, σ, b; evalType = Float64)

Return a list of coefficients for the Laguerre series in Weeks method.

Used in [`GenerateWeeksApproximation`](@ref).
"""
function ComputeCoefficients(lapFunc, N, σ, b; evalType = Float64)

  coefficients = [Calculate_ak(lapFunc, k, N, σ, b, evalType) for k in 0:N]

  return coefficients

end


"""
    Calculate_ak(lapFunc, k, N, σ, b, evalType)

Calculate the `k`th coefficient of the Weeks method approximation using the method by Lyness and Giunta:
Lyness, J. N., & Giunta, G. (1986).
A Modification of the Weeks Method for Numerical Inversion of the Laplace Transform.
Mathematics of Computation, 47(175), 313.
https://doi.org/10.2307/2008097

The underlying integral is performed around a circular contour of radius `r`, which should be less than one, and with `m` quadrature points.

Used in [`ComputeCoefficients`](@ref).
"""
function Calculate_ak(lapFunc, k, N, σ, b, evalType; r = 0.9999, m = 2*N)
  r = 0.9999 #Radius of contour in the method. Should be slightly smaller than 1
  m = 2*N #Number of points in quadrature scheme for the a_k values
  φ = z -> (b/(1-z))*lapFunc( (b/(1-z)) - b/2 + σ  )

  #Calculate a_k
  a_k = real(φ(r)) - imag(φ(r)) + (-1)^k*φ(-r)
  for j in 1:m/2-1
    φ_j = φ( r*exp( 2*convert(evalType,π)*1im*j/m  )  ) #To prevent π being erroneously converted to a Float64, instead of BigFloat, in the relevant case.
    a_k += 2*( real(φ_j)*cos(2*convert(evalType,π)*k*j/m) + imag(φ_j)*sin(2*convert(evalType,π)*k*j/m)   )
  end
  a_k *= (1/(r^k*convert(evalType, m)))

  return a_k

end


end
