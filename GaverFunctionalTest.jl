using ArbNumerics

function FunctionalComparison(func, t, M)

  precisionFactor = 2.1 #See "Multi-precision laplace transform inversion"
  precision = Int(ceil( precisionFactor*M  ))
  if precision < 8
    precision = 8 #ArbNumerics has a minimum decimal precision of 8
  end

  t = ArbFloat(t, digits = precision)
  Dt = typeof(t)




  #||--Gaver functionals--||#
  tau = log(convert(Dt, 2)) / t
  broken = false #used in case of division by zero in wynn rho

  M1 = M #Why?


  G0 = zeros(Dt, M1 + 1)#stores Gaver functionals 0 through M which are all needed for Wynn rho


  G_prev = zeros(Dt, 2*M1 + 1)#stores the previous step in the iteration for the Gaver functionals
  G_next = zeros(Dt, 2*M1 + 1)#stores the next step in the iteration for the Gaver functional


  #Initialise the iteration
  Threads.@threads for n = 1:2M
    samplePoint = ( n*log(convert(Dt,2)) )/(t)
    G_prev[n+1] = samplePoint*func( samplePoint )
  end

  #Calculate Gaver functionals
  for k in 1:M
    Threads.@threads for n in k:2*M-k
      G_next[n+1] = (1+convert(Dt, n//k))*G_prev[n+1] - (convert(Dt, n//k))*G_prev[n+2]
    end
    G_prev = copy(G_next) #Prep for next iteration step
  end

  G0_iter = G_next[1:M1+1]



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

  return [G0_iter  G0]


end


Functionals = FunctionalComparison(Cohen8, 3, 10)
