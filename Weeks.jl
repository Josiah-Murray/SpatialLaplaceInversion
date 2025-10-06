module Weeks

#Define using `WeeksApproximation` function
mutable struct WeeksApproximation
  lapFunc #Laplace domain function to be approximated
  N #Order of the leading term in the polynomial approximation. Should be even.
  σ #Control parameter
  b #Control parameter
  coefficients #Coefficients of the Laguerre series approximating lapFunc
  evalType #The type in which to do the calculations. In general, use Float64. Also built to use BigFloat.
end

#MARK: Main functions

function WeeksApproximation(lapFunc, N, σ, b; evalType = Float64)
  σ = convert(evalType, σ)
  b = convert(evalType, b)
  coefficients = ComputeCoefficients(lapFunc, N, σ, b)
  return WeeksApproximation(lapFunc, N, σ, b, coefficients, evalType)
end

#Uses the clenshaw algorithm to sum the laguerre polynomials
function EvalWeeks(Weeks::WeeksApproximation, t)
  #TODO Implement
  evalType = Weeks.evalType
  N = Weeks.N
  α_n = n -> - (1/(convert(evalType,n)+1))*t + (2*convert(evalType,n) + 1)/(convert(evalType, n)+1)
  β_n = n -> - convert(evalType, n)/(convert(evalType, n)+1)
  p_0 = 1
  p_1 = 1-t

  #Initial terms in the Clenshaw algorithm. Will be updated later to correspond to the two terms before the one we're calculating #TODO: Reword
  bp2 = 0
  bp1 = 0
  b = 0

  for i in N:-1:1
    b = Weeks.coefficients[i+1] + α_n(i)*bp1 + β_n(i+1)*bp2


    bp2 = bp1
    bp1 = b

  end

  approximation = Weeks.coefficients[1] + (1-t)* bp1 + β_n(1)*bp2



  return approximation


end




#MARK: Internal functions


function ComputeCoefficients(lapFunc, N, σ, b; evalType = Float64)

  coefficients = [Calculate_ak(lapFunc, k, N, σ, b, evalType) for k in 0:N]

  return coefficients

end

#Calculates the coefficients of the Weeks method approximation using the method by Lyness and Giunta:
# Lyness, J. N., & Giunta, G. (1986).
# A Modification of the Weeks Method for Numerical Inversion of the Laplace Transform.
# Mathematics of Computation, 47(175), 313.
# https://doi.org/10.2307/2008097
#BUG Not working for k=0 and N, I think
#m: Number of points in quadrature scheme for a_k values
function Calculate_ak(lapFunc, k, N, σ, b, evalType)
  r = 0.9999 #Radius of contour in the method. Should be slightly smaller than 1
  m = 2*N #Number of points in quadrature scheme for the a_k values
  φ = z -> (b/(1-z))*lapFunc( (b/(1-z)) - b/2 + σ  )

  a_k = real(φ(r)) - imag(φ(r)) + (-1)^k*φ(-r)

  #a_k = 0

  for j in 1:m/2-1
    φ_j = φ( r*exp( 2*convert(evalType,π)*1im*j/m  )  ) #To prevent π being erroneously converted to a Float64, instead of BigFloat, in the relevant case.
    a_k += 2*( real(φ_j)*cos(2*convert(evalType,π)*k*j/m) + imag(φ_j)*sin(2*convert(evalType,π)*k*j/m)   )
  end
  a_k *= (1/(r^k*convert(evalType, m)))

  return a_k

end




end
