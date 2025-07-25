#=

This file contains variants of the GWR algorithm.
The purpose is to isolate the effects of changes such as BigFloat vs ArbNumerics
and iterative vs explicit Gaver functional definitions.

=#



#MARK: InverseLaplace.jl
#=

The implementation of GWR in InverseLaplace.jl provided under the license:

"""
Copyright (c) 2015-2018: John Lapeyre.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

Comments are my own.

=#

using ArbNumerics



function gwr_package_GρThread(func, t, M)

  precisionFactor = 2.1 #See "Multi-precision laplace transform inversion"
  precision = Int(ceil( 2.1*M  ))
  if precision < 8
    precision = 8 #ArbNumerics has a minimum decimal precision of 8
  end

  t = ArbFloat(t, digits = precision)
  Dt = typeof(t)




  #||--Gaver functionals--||#
  tau = log(convert(Dt, 2)) / t
  broken = false #used in case of division by zero in wynn rho
  Fi = Array{Dt}(undef, 2 * M)
  @inbounds  for i in 1: 2 * M #The @inbounds improves performance by preventing an internal check.
    Fi[i] = func(i * tau) #Necessary function values are pre-calculated as many are re-used
  end
  M1 = M
  G0 = zeros(Dt, M1 + 1)#stores Gaver functionals 0 through M which are all needed for Wynn rho

  #Calculate Gaver functionals
  Threads.@threads for n in 1:M
    sm = zero(Dt)
    bn = big(n)
    for i in 0:n
      bi = convert(Dt, i)
      sm += convert(Dt, binomial(big(n), big(i))) * (-1)^i * Fi[n + i]
    end
    #The factor bn (k in math) out the front cancels with one of the factorials
    G0[n] = tau * convert(Dt, SpecialFunctions.factorial(2 * bn) /
        (SpecialFunctions.factorial(bn) * SpecialFunctions.factorial(bn - 1))) * sm
  end



    #||--Wynn rho--||#

    #Wynn rho uses values at 2 previous steps. G0 and Gm here.
    Gm = zeros(Dt, M1 + 1)

    #stores current iteration
    Gp = zeros(Dt, M1 + 1)


    best = G0[M1] #Stores approximation at previous step in case of division by zero error.
    for k in 0:M1-2
        Threads.@threads for n in (M1 - 2 - k):-1:0 #Need to go backwards to avoid saving over stuff we need
            expr = G0[n + 2] - G0[n + 1]#Denominator part of recursive step in Wynn rho
            #Check if there will be division by zero. If yes, stop algorithm.
            if expr == 0
                broken = true
                @warn "Division by zero pre-empted on Wynn rho step $k out of $(M1-2).\n    Laplace function: $func"
                break
            end
            #Calculate next terms in Wynn rho
            expr = Gm[n + 2] + (k + 1) / expr
            Gp[n + 1] = expr

            #This true?
            if isodd(k) && n == M1 - 2 - k
                best = expr
            end
        end
        if broken break end
        #Move to next rows of approximations
        for n in 0:(M1-k)
            Gm[n + 1] = G0[n + 1]
            G0[n + 1] = Gp[n + 1]
        end
    end
    best
end




function gwr_package_GThread(func, t, M)

  precisionFactor = 2.1 #See "Multi-precision laplace transform inversion"
  precision = Int(ceil( 2.1*M  ))
  if precision < 8
    precision = 8 #ArbNumerics has a minimum decimal precision of 8
  end

  t = ArbFloat(t, digits = precision)
  Dt = typeof(t)




  #||--Gaver functionals--||#
  tau = log(convert(Dt, 2)) / t
  broken = false #used in case of division by zero in wynn rho
  Fi = Array{Dt}(undef, 2 * M)
  @inbounds  for i in 1: 2 * M #The @inbounds improves performance by preventing an internal check.
    Fi[i] = func(i * tau) #Necessary function values are pre-calculated as many are re-used
  end
  M1 = M
  G0 = zeros(Dt, M1 + 1)#stores Gaver functionals 0 through M which are all needed for Wynn rho

  #Calculate Gaver functionals
  Threads.@threads for n in 1:M
    sm = zero(Dt)
    bn = big(n)
    for i in 0:n
      bi = convert(Dt, i)
      sm += convert(Dt, binomial(big(n), big(i))) * (-1)^i * Fi[n + i]
    end
    #The factor bn (k in math) out the front cancels with one of the factorials
    G0[n] = tau * convert(Dt, SpecialFunctions.factorial(2 * bn) /
        (SpecialFunctions.factorial(bn) * SpecialFunctions.factorial(bn - 1))) * sm
  end



    #||--Wynn rho--||#

    #Wynn rho uses values at 2 previous steps
    Gm = zeros(Dt, M1 + 1)
    Gp = zeros(Dt, M1 + 1)


    best = G0[M1] #Stores approximation at previous step in case of division by zero error.
    for k in 0:M1-2
        for n in (M1 - 2 - k):-1:0 #Need to go backwards to avoid saving over stuff we need
            expr = G0[n + 2] - G0[n + 1]#Denominator part of recursive step in Wynn rho
            #Check if there will be division by zero. If yes, stop algorithm.
            if expr == 0
                broken = true
                @warn "Division by zero pre-empted on Wynn rho step $k out of $(M1-2).\n    Laplace function: $func"
                break
            end
            #Calculate next terms in Wynn rho
            expr = Gm[n + 2] + (k + 1) / expr
            Gp[n + 1] = expr

            #This true?
            if isodd(k) && n == M1 - 2 - k
                best = expr
            end
        end
        if broken break end
        #Move to next rows of approximations
        for n in 0:(M1-k)
            Gm[n + 1] = G0[n + 1]
            G0[n + 1] = Gp[n + 1]
        end
    end
    best
end




function gwr_package_arbnumerics(func, t, M)

  precisionFactor = 2.1 #See "Multi-precision laplace transform inversion"
  precision = Int(ceil( 2.1*M  ))
  if precision < 8
    precision = 8 #ArbNumerics has a minimum decimal precision of 8
  end

  t = ArbFloat(t, digits = precision)
  Dt = typeof(t)




  #||--Gaver functionals--||#
  tau = log(convert(Dt, 2)) / t
  broken = false #used in case of division by zero in wynn rho
  Fi = Array{Dt}(undef, 2 * M)
  @inbounds  for i in 1: 2 * M #The @inbounds improves performance by preventing an internal check.
    Fi[i] = func(i * tau) #Necessary function values are pre-calculated as many are re-used
  end
  M1 = M
  G0 = zeros(Dt, M1 + 1)#stores Gaver functionals 0 through M which are all needed for Wynn rho

  #Calculate Gaver functionals
  for n in 1:M
    sm = zero(Dt)
    bn = big(n)
    for i in 0:n
      bi = convert(Dt, i)
      sm += convert(Dt, binomial(big(n), big(i))) * (-1)^i * Fi[n + i]
    end
    #The factor bn (k in math) out the front cancels with one of the factorials
    G0[n] = tau * convert(Dt, SpecialFunctions.factorial(2 * bn) /
        (SpecialFunctions.factorial(bn) * SpecialFunctions.factorial(bn - 1))) * sm
  end



    #||--Wynn rho--||#

    #Wynn rho uses values at 2 previous steps
    Gm = zeros(Dt, M1 + 1)
    Gp = zeros(Dt, M1 + 1)


    best = G0[M1] #Stores approximation at previous step in case of division by zero error.
    for k in 0:M1-2
        for n in (M1 - 2 - k):-1:0 #Need to go backwards to avoid saving over stuff we need
            expr = G0[n + 2] - G0[n + 1]#Denominator part of recursive step in Wynn rho
            #Check if there will be division by zero. If yes, stop algorithm.
            if expr == 0
                broken = true
                @warn "Division by zero pre-empted on Wynn rho step $k out of $(M1-2).\n    Laplace function: $func"
                break
            end
            #Calculate next terms in Wynn rho
            expr = Gm[n + 2] + (k + 1) / expr
            Gp[n + 1] = expr

            #This true?
            if isodd(k) && n == M1 - 2 - k
                best = expr
            end
        end
        if broken break end
        #Move to next rows of approximations
        for n in 0:(M1-k)
            Gm[n + 1] = G0[n + 1]
            G0[n + 1] = Gp[n + 1]
        end
    end
    best
end



function gwr_package(func, t, M)
    Dt = typeof(t)
    bM = convert(Dt, M) #Unused




    #||--Gaver functionals--||#
    tau = log(convert(Dt, 2)) / t
    broken = false #used in case of division by zero in wynn rho
    Fi = Array{Dt}(undef, 2 * M)
@inbounds  for i in 1: 2 * M #The @inbounds improves performance by preventing an internal check.
        Fi[i] = func(i * tau) #Necessary function values are pre-calculated as many are re-used
    end
    M1 = M
    G0 = zeros(Dt, M1 + 1)#stores Gaver functionals 0 through M which are all needed for Wynn rho

    #Calculate Gaver functionals
    for n in 1:M
        sm = zero(Dt)
        bn = big(n)
        for i in 0:n
            bi = convert(Dt, i)
            sm += binomial(big(n), big(i)) * (-1)^i * Fi[n + i]
        end
        #The factor bn (k in math) out the front cancels with one of the factorials
        G0[n] = tau * SpecialFunctions.factorial(2 * bn) /
            (SpecialFunctions.factorial(bn) * SpecialFunctions.factorial(bn - 1)) * sm
    end



    #||--Wynn rho--||#

    #Wynn rho uses values at 2 previous steps
    Gm = zeros(Dt, M1 + 1)
    Gp = zeros(Dt, M1 + 1)


    best = G0[M1] #Stores approximation at previous step in case of division by zero error.
    for k in 0:M1-2
        for n in (M1 - 2 - k):-1:0 #Need to go backwards to avoid saving over stuff we need
            expr = G0[n + 2] - G0[n + 1]#Denominator part of recursive step in Wynn rho
            #Check if there will be division by zero. If yes, stop algorithm.
            if expr == 0
                broken = true
                @warn "Division by zero pre-empted on Wynn rho step $k out of $(M1-2).\n    Laplace function: $func"
                break
            end
            #Calculate next terms in Wynn rho
            expr = Gm[n + 2] + (k + 1) / expr
            Gp[n + 1] = expr

            #This true?
            if isodd(k) && n == M1 - 2 - k
                best = expr
            end
        end
        if broken break end
        #Move to next rows of approximations
        for n in 0:(M1-k)
            Gm[n + 1] = G0[n + 1]
            G0[n + 1] = Gp[n + 1]
        end
    end
    best
end
