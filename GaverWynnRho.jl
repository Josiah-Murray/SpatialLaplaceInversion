#This file contains the functions for my implementation of the Gaver-Wynn rho algrothim.

#The arbitrary precision is handled by the ArbNumerics package.
#Notably, because the Gaver functionals involve the evaluation of the Laplace domain function,
# it should also be designed to use arbitrary precision in its calculations.




#The implementation is based on
#=
Abate, J., & Valkó, P. P. (2004). Multi‐precision Laplace transform inversion.
International Journal for Numerical Methods in Engineering, 60(5), 979–993.
https://doi.org/10.1002/nme.995
=#

#The GWR method requires careful handling of precision in the calculations
#to avoid cancellation errors.
using ArbNumerics


#=
The primary function which runs the inversion.
Inverts a function 'LaplaceFunction' of the form F(s), for some s, and returns
f(t) for the given 't'. M is the number of terms in the approximation.
See the referrenced paper at the top of this file for more information.
=#
function gwr( LaplaceFunction, t, M  )


  #As per Abate and Valkó (2004), the precision should be set to (2.1)M
  precision = Int(ceil(2.1*M))
  t = ArbReal(t, digits = precision)

  #Calculate array of Gaver functionals
  f_k_vector = [GaverFunctional(LaplaceFunction, t, k, precision) for k in 0:M ]


  #Apply Wynn rho algorithm to approximate inversion.
  ft = WynnRho(f_k_vector, precision)

  return ft

end#function




#|||| DEPRECATED||||#
#Reliance on 'binomial()' causes problems
#Calculates the kth Gaver functional for a Laplace-domain function 'LaplaceFunction
#a time 't' and with precision.
function GaverFunctional(LaplaceFunction, t::ArbReal, k::Int, precision::Int)

  #Alias the typing to avoid excessively long lines of code.
  A(n) = ArbReal(n, digits = precision)


  #BUG binomial can only handle values up to a certain size.
  #Calculate the functional f_k based on the formula (See equation (4) in referenced paper).
  f_k = ArbReal(0, digits = precision)
  for j = 0:k
    f_k += (-1)^j*binomial(k,j)*LaplaceFunction( (k+j)*log(A(2))/t  )
  end
  f_k *= (k*log(A(2))/t)*binomial(2k,k)

  #A one-line version
  #f_k = (k*log(A(2))/t)*binomial(2k,k)*sum( [(-1)^j*binomial(k,j)*LaplaceFunction( (k+j)*log(A(2))/t  ) for j in 0:k ] )

  return f_k

end


#Applies the Wynn rho algorithm to approximate the solution to the inversion.
#Takes an vector of Functionals from 'GaverFunctional' and a precision (integer),
# and outputs the approximation to the inversion.
function WynnRho(FunctionalVector, precision)
  M = length(FunctionalVector)-1

  #stores the elements in the Wynn rho sequence.
  #In particular, ρ_k^(n) is stored at position ρ[n+1,k+2].
  #Note that in the definition of the algorithm, k begins at -1 and
  #n begins at 0.
  #To make it easier to read, I have left the +1 and +2 unsimplified.
  #That means ρ[0+1, 2+2], for instance, is equivalent to ρ_2^(0).
  ρ = fill(ArbReal(0, digits = precision), M+1, M+2)

  #First column only zeros, then second column is populated by the functionals.
  ρ[:, 2] = FunctionalVector

  #Do iteration based on formula (See equation (6) in referenced paper.)
  for k = 1:M
    for n = 0:(M-k)
      ρ[n+1,k+2] = ρ[n+1+1, k-2+2] + k/(ρ[n+1+1, k-1+2] - ρ[n+1, k-1+2])
    end
  end

  #Return element corresponding to the solution (See equation (7) in referenced paper.)
  return ρ[0+1,M+2]


end
