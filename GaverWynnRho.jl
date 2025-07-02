#This file contains the functions for my implementation of the Gaver-Wynn rho algrothim.

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
  precision = ceil(2.1*M)

  #Calculate array of Gaver functionals

  #Apply Wynn rho algorithm





end#function




function GaverFunctional(LaplaceFunction, t::ArbReal, k:::Int, precision::Int)

  #Alias the typing to avoid excessively long lines of code.
  A(n) = ArbReal(n, digits = precision)

  f_k = ArbReal(0, digits = precision)
  for j = 1:k
    #TODO Double check for fatigue errors
    f_k += (-1)^j*binomial(k,j)*LaplaceFunction( (k+j)*log(A(2))/t  )
  end

  f_k *= k*log(A(2))/t

  return f_k

end
