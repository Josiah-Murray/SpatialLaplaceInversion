#This file contains a set of test functions for evaluating the performance of
#numerical Laplace inversion schemes


using SpecialFunctions

#||||---Cohen----||||#

#Cohen, A. M. (2007).
#Numerical Methods for Laplace Transform Inversion (1st ed., Vol. 5).
#Springer US.
#https://doi.org/10.1007/978-0-387-68855-8

function Cohen1(s)
  return (s^2 + 1)^(-1/2)
end

function Cohen1_exact(t)
  return besselj0(t)
end

function Cohen2(s)
  return s^(-1/2)*exp(-1/s)
end

function Cohen2_exact(t)
  return (π*t)^(-1/2)*cos(2*sqrt(t))
end

function Cohen3(s)
  return (s+1/2)^(-1)
end

function Cohen3_exact(t)
  return exp(-t/2)
end

function Cohen4(s)
  return 1/( (s+0.2)^2 + 1)
end

function Cohen4_exact(s)
  return exp(-0.2t)*sin(t)
end

function Cohen5(s)
  return s^(-1)
end

function Cohen5_exact(t)
  return 1
end

function Cohen6(s)
  return s^(-2)
end

function Cohen6_exact(t)
  return t
end

function Cohen7(s)
  return (s+1)^(-2)
end

function Cohen7_exact(t)
  return t*exp(-t)
end

function Cohen8(s)
  return (s^2 + 1)^(-1)
end


function Cohen8_exact(t)
  return sin(t)
end

function Cohen9(s)
  return s^(-1/2)
end

function Cohen9_exact(t)
  return (πt)^(-1/2)
end

function Cohen10(s)
  return s^(-1)*exp(-5s)
end

function Cohen10_exact(t)
  if t<5
    return 0
  elseif t>5
    return 1
  else
    return 0.5
  end
end

function Cohen11(s)
  return s^(-1)*log(s)
end

function Cohen11_exact(t)
  return -Base.MathConstants.eulergamma - log(t)
end

function Cohen12(s)
  return (s(1+exp(-s)))^(-1)
end

function Cohen12_exact(t)
  if isinteger(t)
    return 0.5
  elseif mod(floor(t),2) == 0
    return 1
  else
    return 0
  end
end


function Cohen13(s)
  return (s^2-1)(s^2+1)^(-2)
end

function Cohen13_exact(t)
  return t*cos(t)
end

function Cohen14(s)
  return (s+1/2)^(1/2) - (s+1/4)^(1/2)
end

function Cohen14_exact(t)
  return (exp(-t/3) - exp(-t/2))*(4*π*t^3)^(-1/2)
end


function Cohen15(s)
  return exp(-4*s^(1/2))
end

function Cohen15_exact(t)
  return 2*exp(-4/t)*(π*t^3)^(-1/2)
end

function Cohen16(s)
  return atan(1/s)
end

function Cohen16_exact(t)
  return t^(-1)*sin(t)
end

function Cohen17(s)
  return 1/(s^3)
end

function Cohen17_exact(t)
  return (1/2)*t^2
end

function Cohen18(s)
  return 1/(s^2 + s + 1)
end

function Cohen18_exact(t)
  return (2/sqrt(3))*exp(-t/2)*sin(sqrt(3)*t/2)
end

function Cohen19(s)
  return 3/(s^2 - 9)
end

function Cohen19_exact(t)
  return sinh(3t)
end

function Cohen20(s)
  return 120/(s^6)
end

function Cohen20_exact(t)
  return t^5
end

function Cohen21(s)
  return s/((s^2 + 1)^2)
end

function Cohen21_exact(t)
  return (1/2)*t*sin(t)
end

function Cohen22(s)
  return (s+1)^(-1) - (s+1000)^(-1)
end

function Cohen22_exact(t)
  return exp(-t) + exp(-1000t)
end

function Cohen23(s)
  return s/(s^2+1)
end

function Cohen23_exact(t)
  return cos(t)
end

function Cohen24(s)
  return 1/((s-0.25)^2)
end

function Cohen24_exact(t)
  return t*exp(t/4)
end

function Cohen25(s)
  return 1/(s*sqrt(s))
end

function Cohen25_exact(t)
  return 2*sqrt(t/π)
end

function Cohen26(s)
  return 1/((s+1)^(1/2))
end

function Cohen26_exact(t)
  return exp(-t)/sqrt(π*t)
end

function Cohen27(s)
  return (s+2)/(s*sqrt(s))
end

function Cohen27_exact(t)
  return  (1+4*t)/sqrt(π*t)
end


function Cohen28(s)
  return 1/((s^2+1)^2)
end

function Cohen28_exact(t)
  return (1/2)*(sin(t) - t*cos(t))
end

function Cohen29(s)
  return 1/(s*(s+1)^2)
end

function Cohen29_exact(t)
  return 1-exp(-t)*(1+t)
end

function Cohen30(s)
  return 1/(s^3 - 8)
end

function Cohen30_exact(t)
  return (1/12)*exp(-t)*( exp(3t) - cos(sqrt(3)*t) - sqrt(3)sin(sqrt(3)*t)  )
end

function Cohen31(s)
  return log( (s^2+1)/(s^2+4)  )
end

function Cohen31_exact(t)
  return 2(cos(2t) - cos(t))/t
end

function Cohen32(s)
  return log((s+1)/s)
end

function Cohen32_exact(t)
  return (1-exp(-t))/t
end

function Cohen33(s)
  return (1-exp(-s))/(s^2)
end

function Cohen33_exact(t)
  if 0 <= t <= 1
    return t
  else
    return 0
  end
end

function Cohen34(s)
  return 1/(s*(1+exp(s)))
end

function Cohen34_exact(t)
  if isinteger(t)
    return 0.5
  end
  if mod(floor(t), 2) == 0
    return 0
  else
    return 1
  end
end


function Cohen35(s)
  return 1/(s^(1/2) + s^(1/3))
end

#The exact inversion of Cohen35 is referenced as appearing
# in Section 4.1
#It actually appears as Example 3.6
#=
N=1000 is chosen somewhat arbitrarily. Numerical experiments indicate that
for t≈10, the terms are approximately of the order 10^-132 for n≈1000. which
is far beyond the accuracies usually recorded for numerical Laplace inversion
algorithms.
=#
function Cohen35_exact(t)
  N = 600
  return t^(1/2)*sum(  [ (  (-1)^n*t^(n/6)/gamma((n+3)/6)  )   for n ∈ 0:N ]   )
end
