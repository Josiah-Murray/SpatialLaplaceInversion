module GaverWynnRho

using SpecialFunctions


#TODO: Add comment and citation to InverseLaplace.jl
#TODO: Is it compatible with the ArbNumerics?
#TODO: Implement shifting.
function GWR(func, t, M)

  Dt = typeof(t)
  Dt = Dt <: Complex ? Dt : Complex{Dt}
  if Dt <: Int
      @error "Type conversions fail for integers. Consider using a t with type BigFloat."
    end

  #||--Gaver functionals--||#
  tau = log(convert(Dt, 2)) / t #Sample points of Laplace domain function, `func` for Gaver functionals.
  Fi = Array{Dt}(undef, 2 * M) #Stores evaluations of `func`
  @inbounds  for i in 1: 2 * M #The @inbounds improves performance by preventing an internal check. Will over right other data if implemented wrong.
    Fi[i] = func(i * tau) #Necessary function values are pre-calculated as many are re-used
  end

  ρ_0 = zeros(Dt, M + 1)#stores Gaver functionals 0 through M which are all needed for Wynn rho. Later over-written as part of the Wynn rho algorithm.

  #Calculate Gaver functionals
  for n in 1:M
    sm = zero(Dt)
    bn = big(n)
    for i in 0:n
      sm += binomial(bn, big(i)) * (-1)^i * Fi[n + i]
    end
    #The factor bn (k in math) out the front cancels with one of the factorials
    ρ_0[n+1] = tau * SpecialFunctions.factorial(2 * bn) /
      (SpecialFunctions.factorial(bn) * SpecialFunctions.factorial(bn - 1)) * sm
  end


  #||--Wynn rho--||#


  broken = false #used in case of division by zero in wynn rho
  #Wynn rho uses values at 2 previous steps
  ρ_m = zeros(Dt, M + 1)#Stores the ρ values for ρ_{i-1}^( ... )
  ρ_p = zeros(Dt, M + 1)#Stores the newly calculated ρ values for ρ_{i+1}^( ... )
  #ρ_{i}^( ... ) are stored in `ρ_0`


  best_approximation = ρ_0[M] #Stores last valid approximation in case of division by zero. This can happen when the method recovers the exact solution.
  for k in 1:M

    for n in (M - k):-1:0
      denominator = ρ_0[n + 2] - ρ_0[n + 1]#Denominator part of recursive step in Wynn rho. Check if zero before proceeding.
      #Check if there will be division by zero. If yes, stop algorithm.
      if denominator == 0
        broken = true
        @warn "Division by zero pre-empted on Wynn rho step $k out of $(M-2).\n    Laplace function: $func"
        break
      end
      #Calculate next terms in Wynn rho
      ρ_p[n + 1] = ρ_m[n + 2] + (k) / denominator

      #Check if valid approximation
      if iseven(k) && n == M - k
        best = denominator
      end
    end


    if broken
      break
    end

    #Move to next rows of approximations
    for n in 0:(M-k)
      ρ_m[n + 1] = ρ_0[n + 1]
      ρ_0[n + 1] = ρ_p[n + 1]
    end
  end
  best_approximation
end

#TODO: Consider moving to a different file.
#TODO: Fix naming of variables
#func_array(s) should return an array of values (possible multiple dimensions).
function GWR_array(input_func_array, t, M; shift_parameter = 0)
    Dt = typeof(t)
    Dt = Dt <: Complex ? Dt : Complex{Dt}
    if Dt <: Int
      @error "Type conversions fail for integers. Consider using a t with type BigFloat."
    end
    bM = convert(Dt, M) #Unused
    #println("Update")#BUG Debug line

    func_array = s -> input_func_array(s+shift_parameter)


    #||--Gaver functionals--||#
    tau = log(convert(Dt, 2)) / t
    broken = false #used in case of division by zero in wynn rho
    Fi = Array{Array{Dt}}(undef, 2 * M)
  @inbounds  for i in 1: 2 * M #The @inbounds improves performance by preventing an internal check. Will over right other data if implemented wrong.
        Fi[i] = func_array(i * tau) #Necessary function values are pre-calculated as many are re-used
    end

    M1 = M
    #G0 = Array{typeof(Array{Dt}(undef, size(Fi[1])...))}(undef, M1+1) #stores Gaver functionals 0 through M which are all needed for Wynn rho
    G0 = [zeros(Dt, size(Fi[1])...) for i in 1:(M1+1)]#fi[1] is a stand in for the dimensions of func_arrray.
    #Calculate Gaver functionals
    for n in 1:M
        #sm = Array{Dt}(undef, size(Fi[1])...) #Assumes all outputs of func_array are the same size.
        sm = zeros(Dt, size(Fi[1])...)
        bn = big(n)
        for i in 0:n
            bi = convert(Dt, i)
            sm += binomial(big(n), big(i)) * (-1)^i * Fi[n + i]
        end
        #The factor bn (k in math) out the front cancels with one of the factorials
        G0[n+1] = tau * SpecialFunctions.factorial(2 * bn) /
            (SpecialFunctions.factorial(bn) * SpecialFunctions.factorial(bn - 1)) * sm
    end




    #||--Wynn rho--||#

    #Wynn rho uses values at 2 previous steps
    ArrayType = typeof(Array{Dt}(undef, size(Fi[1])...))
    Gm = [zeros(Dt, size(Fi[1])...) for i in 1:(M1+1)]#fi[1] is a stand in for the dimensions of func_arrray.
    Gp = [zeros(Dt, size(Fi[1])...) for i in 1:(M1+1)]

    mask = trues(size(Fi[1])...) #Used to 'mask out' any of the values where we have stopped the calculation early.



    best = G0[M1] #Stores approximation at previous step in case of division by zero error.
    for k in 1:M1
        for n in (M1 - k):-1:0
            expr = G0[n + 2] - G0[n + 1]#Denominator part of recursive step in Wynn rho
            #Check if there will be division by zero. If yes, stop algorithm.
            if 0 ∈ expr #BUG Could this be done more efficiently?
                mask_indices = findall(x -> x==0, expr)

                #BUG: Implement mask update
                for indices in mask_indices
                  mask[indices] = false
                end
                #@warn "Division by zero pre-empted on Wynn rho step $k out of $(M1-2).\n    Laplace function: $func_array \n Mask: $mask"
                if true ∉ mask #If all have been masked out, stop the procedure.
                  break
                end

            end

            #Calculate next terms in Wynn rho
            expr = Gm[n + 2] + (k) ./ expr
            Gp[n + 1] = expr #TODO: Implement masking

            #This true?
            if iseven(k) && n == M1 - k
              for i in eachindex(expr)
                if(mask[i])
                  best[i] = expr[i]
                end
              end
            end
        end

        #Move to next rows of approximations
        for n in 0:(M1-k)
            Gm[n + 1] = G0[n + 1] #TODO Implement masking
            G0[n + 1] = Gp[n + 1]

            #Gm = copy(G0)
            #G0 = copy(Gp)

        end
    end


    return best*exp(shift_parameter)
end





end
