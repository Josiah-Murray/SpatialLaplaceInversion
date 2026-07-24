module TestFunctions

#This file contains a set of test functions for evaluating the performance of
#numerical Laplace inversion schemes


using ..SpecialFunctions
#||||---Cohen----||||#

#Cohen, A. M. (2007).
#Numerical Methods for Laplace Transform Inversion (1st ed., Vol. 5).
#Springer US.
#https://doi.org/10.1007/978-0-387-68855-8


#Takes in a function, 'inversionScheme', for inverting the Laplace transform and
#returns a 2 by 35 array where the rows correspond to the tests (i.e. row
#one corresponds to 'Cohen1' and 'Cohen1_exact),
#the first column corresponds to the exact solutions, and the second to the
#approximate solutions.
function CohenSuite( inversionScheme, t)
    exact = Vector{Float64}(undef, 35)
    approx = Vector{Float64}(undef, 35)

    exact[1] = Cohen1_exact(t)
    approx[1] = inversionScheme(Cohen1, t)

    exact[2] = Cohen2_exact(t)
    approx[2] = inversionScheme(Cohen2, t)

    exact[3] = Cohen3_exact(t)
    approx[3] = inversionScheme(Cohen3, t)

    exact[4] = Cohen4_exact(t)
    approx[4] = inversionScheme(Cohen4, t)

    exact[5] = Cohen5_exact(t)
    approx[5] = inversionScheme(Cohen5, t)

    exact[6] = Cohen6_exact(t)
    approx[6] = inversionScheme(Cohen6, t)

    exact[7] = Cohen7_exact(t)
    approx[7] = inversionScheme(Cohen7, t)

    exact[8] = Cohen8_exact(t)
    approx[8] = inversionScheme(Cohen8, t)

    exact[9] = Cohen9_exact(t)
    approx[9] = inversionScheme(Cohen9, t)

    exact[10] = Cohen10_exact(t)
    approx[10] = inversionScheme(Cohen10, t)

    exact[11] = Cohen11_exact(t)
    approx[11] = inversionScheme(Cohen11, t)

    exact[12] = Cohen12_exact(t)
    approx[12] = inversionScheme(Cohen12, t)

    exact[13] = Cohen13_exact(t)
    approx[13] = inversionScheme(Cohen13, t)

    exact[14] = Cohen14_exact(t)
    approx[14] = inversionScheme(Cohen14, t)

    exact[15] = Cohen15_exact(t)
    approx[15] = inversionScheme(Cohen15, t)

    exact[16] = Cohen16_exact(t)
    approx[16] = inversionScheme(Cohen16, t)

    exact[17] = Cohen17_exact(t)
    approx[17] = inversionScheme(Cohen17, t)

    exact[18] = Cohen18_exact(t)
    approx[18] = inversionScheme(Cohen18, t)

    exact[19] = Cohen19_exact(t)
    approx[19] = inversionScheme(Cohen19, t)

    exact[20] = Cohen20_exact(t)
    approx[20] = inversionScheme(Cohen20, t)

    exact[21] = Cohen21_exact(t)
    approx[21] = inversionScheme(Cohen21, t)

    exact[22] = Cohen22_exact(t)
    approx[22] = inversionScheme(Cohen22, t)

    exact[23] = Cohen23_exact(t)
    approx[23] = inversionScheme(Cohen23, t)

    exact[24] = Cohen24_exact(t)
    approx[24] = inversionScheme(Cohen24, t)

    exact[25] = Cohen25_exact(t)
    approx[25] = inversionScheme(Cohen25, t)

    exact[26] = Cohen26_exact(t)
    approx[26] = inversionScheme(Cohen26, t)

    exact[27] = Cohen27_exact(t)
    approx[27] = inversionScheme(Cohen27, t)

    exact[28] = Cohen28_exact(t)
    approx[28] = inversionScheme(Cohen28, t)

    exact[29] = Cohen29_exact(t)
    approx[29] = inversionScheme(Cohen29, t)

    exact[30] = Cohen30_exact(t)
    approx[30] = inversionScheme(Cohen30, t)

    exact[31] = Cohen31_exact(t)
    approx[31] = inversionScheme(Cohen31, t)

    exact[32] = Cohen32_exact(t)
    approx[32] = inversionScheme(Cohen32, t)

    exact[33] = Cohen33_exact(t)
    approx[33] = inversionScheme(Cohen33, t)

    exact[34] = Cohen34_exact(t)
    approx[34] = inversionScheme(Cohen34, t)

    exact[35] = Cohen35_exact(t)
    approx[35] = inversionScheme(Cohen35, t)

    return [exact approx]
end


#Same as 'CohenSuite' except it excludes solution functions which aren't continuously differentiable
#and solution functions which have undefined limits as t approaches 0. Namely, 2,9, 10, 11, 12, 14, 15, 26, 33, and 34
function CohenSuiteContDiff( inversionScheme, t)
  exact = Vector{Float64}(undef, 25)
    approx = Vector{Float64}(undef, 25)
    function_array = Vector{Function}(undef, 25)
    #TODO add way to see which functions are used.

    i=1
    exact[i] = Cohen1_exact(t)
    approx[i] = inversionScheme(Cohen1, t)
    function_array[i] = Cohen1

    #=
    exact[2] = Cohen2_exact(t)
    approx[2] = inversionScheme(Cohen2, t)
    =#
    i+=1
    exact[i] = Cohen3_exact(t)
    approx[i] = inversionScheme(Cohen3, t)
    function_array[i] = Cohen3

    i+=1
    exact[i] = Cohen4_exact(t)
    approx[i] = inversionScheme(Cohen4, t)
    function_array[i] = Cohen4


    i+=1
    exact[i] = Cohen5_exact(t)
    approx[i] = inversionScheme(Cohen5, t)
    function_array[i] = Cohen5

    i+=1
    exact[i] = Cohen6_exact(t)
    approx[i] = inversionScheme(Cohen6, t)
    function_array[i] = Cohen6

    i+=1
    exact[i] = Cohen7_exact(t)
    approx[i] = inversionScheme(Cohen7, t)
    function_array[i] = Cohen7

    i+=1
    exact[i] = Cohen8_exact(t)
    approx[i] = inversionScheme(Cohen8, t)
    function_array[i] = Cohen8

    #=
    exact[9] = Cohen9_exact(t)
    approx[9] = inversionScheme(Cohen9, t)
    =#


    #=
    exact[10] = Cohen10_exact(t)
    approx[10] = inversionScheme(Cohen10, t)
    =#

    #=
    exact[11] = Cohen11_exact(t)
    approx[11] = inversionScheme(Cohen11, t)
    =#

    #=
    exact[12] = Cohen12_exact(t)
    approx[12] = inversionScheme(Cohen12, t)
    =#

    i+=1
    exact[i] = Cohen13_exact(t)
    approx[i] = inversionScheme(Cohen13, t)
    function_array[i] = Cohen13

    #=
    exact[14] = Cohen14_exact(t)
    approx[14] = inversionScheme(Cohen14, t)
    =#

    #=
    exact[15] = Cohen15_exact(t)
    approx[15] = inversionScheme(Cohen15, t)
    =#

    i+=1
    exact[i] = Cohen16_exact(t)
    approx[i] = inversionScheme(Cohen16, t)
    function_array[i] = Cohen16

    i+=1
    exact[i] = Cohen17_exact(t)
    approx[i] = inversionScheme(Cohen17, t)
    function_array[i] = Cohen17

    i+=1
    exact[i] = Cohen18_exact(t)
    approx[i] = inversionScheme(Cohen18, t)
    function_array[i] = Cohen18

    i+=1
    exact[i] = Cohen19_exact(t)
    approx[i] = inversionScheme(Cohen19, t)
    function_array[i] = Cohen19

    i+=1
    exact[i] = Cohen20_exact(t)
    approx[i] = inversionScheme(Cohen20, t)
    function_array[i] = Cohen20

    i+=1
    exact[i] = Cohen21_exact(t)
    approx[i] = inversionScheme(Cohen21, t)
    function_array[i] = Cohen21

    i+=1
    exact[i] = Cohen22_exact(t)
    approx[i] = inversionScheme(Cohen22, t)
    function_array[i] = Cohen22

    i+=1
    exact[i] = Cohen23_exact(t)
    approx[i] = inversionScheme(Cohen23, t)
    function_array[i] = Cohen23

    i+=1
    exact[i] = Cohen24_exact(t)
    approx[i] = inversionScheme(Cohen24, t)
    function_array[i] = Cohen24

    i+=1
    exact[i] = Cohen25_exact(t)
    approx[i] = inversionScheme(Cohen25, t)
    function_array[i] = Cohen25

    #=
    exact[26] = Cohen26_exact(t)
    approx[26] = inversionScheme(Cohen26, t)
    =#

    i+=1
    exact[i] = Cohen27_exact(t)
    approx[i] = inversionScheme(Cohen27, t)
    function_array[i] = Cohen27

    i+=1
    exact[i] = Cohen28_exact(t)
    approx[i] = inversionScheme(Cohen28, t)
    function_array[i] = Cohen28

    i+=1
    exact[i] = Cohen29_exact(t)
    approx[i] = inversionScheme(Cohen29, t)
    function_array[i] = Cohen29

    i+=1
    exact[i] = Cohen30_exact(t)
    approx[i] = inversionScheme(Cohen30, t)
    function_array[i] = Cohen30

    i+=1
    exact[i] = Cohen31_exact(t)
    approx[i] = inversionScheme(Cohen31, t)
    function_array[i] = Cohen31

    i+=1
    exact[i] = Cohen32_exact(t)
    approx[i] = inversionScheme(Cohen32, t)
    function_array[i] = Cohen32

    #=
    exact[33] = Cohen33_exact(t)
    approx[33] = inversionScheme(Cohen33, t)
    =#

    #=
    exact[34] = Cohen34_exact(t)
    approx[34] = inversionScheme(Cohen34, t)
    =#

    i+=1
    exact[i] = Cohen35_exact(t)
    approx[i] = inversionScheme(Cohen35, t)
    function_array[i] = Cohen35

    return [exact approx function_array]
end



#----#


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

function Cohen4_exact(t)
  return exp(-0.2*t)*sin(t)
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
  return (π*t)^(-1/2)
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
  return (s*(1+exp(-s)))^(-1)
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
  return (s^2-1)*(s^2+1)^(-2)
end

function Cohen13_exact(t)
  return t*cos(t)
end

function Cohen14(s)
  return (s+1/2)^(1/2) - (s+1/4)^(1/2)
end

function Cohen14_exact(t)
  return (exp(-t/4) - exp(-t/2))*(4*π*t^3)^(-1/2)
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
    return 1
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
I've chosen N=1000 chosen somewhat arbitrarily. Numerical experiments indicate that
for t≈10, the terms are approximately of the order 10^-132 for n≈1000. which
is far beyond the accuracies usually recorded for numerical Laplace inversion
algorithms.
=#
function Cohen35_exact(t)
  N = 600
  return t^(-1/2)*sum(  [ (  (-1)^n*t^(n/6)/gamma((n+3)/6)  )   for n ∈ 0:N ]   )
end


end #module
