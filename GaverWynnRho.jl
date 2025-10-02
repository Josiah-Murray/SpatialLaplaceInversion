module GaverWynnRho

using SpecialFunctions


#TODO: Rewrite to make more consistent with standard notation.
function GWR(func, t, M)
    Dt = typeof(t)
    bM = convert(Dt, M) #Unused




    #||--Gaver functionals--||#
    tau = log(convert(Dt, 2)) / t
    broken = false #used in case of division by zero in wynn rho
    Fi = Array{Dt}(undef, 2 * M)
  @inbounds  for i in 1: 2 * M #The @inbounds improves performance by preventing an internal check. Will over right other data if implemented wrong.
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
        G0[n+1] = tau * SpecialFunctions.factorial(2 * bn) /
            (SpecialFunctions.factorial(bn) * SpecialFunctions.factorial(bn - 1)) * sm
    end


    #BUG: G0[end] is always equal to zero for some reason.
    #I think G0[n] should actually be G0[n+1].
    ##TODO: Don't forget to take out the @inbounds while testing
    #print(G0[end])



    #||--Wynn rho--||#

    #Wynn rho uses values at 2 previous steps
    Gm = zeros(Dt, M1 + 1)
    Gp = zeros(Dt, M1 + 1)


    best = G0[M1] #Stores approximation at previous step in case of division by zero error.
    for k in 1:M1
        for n in (M1 - k):-1:0
            expr = G0[n + 2] - G0[n + 1]#Denominator part of recursive step in Wynn rho
            #Check if there will be division by zero. If yes, stop algorithm.
            if expr == 0
                broken = true
                @warn "Division by zero pre-empted on Wynn rho step $k out of $(M1-2).\n    Laplace function: $func"
                break
            end
            #Calculate next terms in Wynn rho
            expr = Gm[n + 2] + (k) / expr
            Gp[n + 1] = expr

            #This true?
            if iseven(k) && n == M1 - k
                best = expr
            end
        end
        if broken break end
        #Move to next rows of approximations
        for n in 0:(M1-k)
            Gm[n + 1] = G0[n + 1]
            G0[n + 1] = Gp[n + 1]

            #Gm = copy(G0)
            #G0 = copy(Gp)

        end
    end
    best
end



#TODO: Implement.
#func_array(s) should return an array of values (possible multiple dimensions).
function GWR_array(func_array, t, M)
    Dt = typeof(t)
    bM = convert(Dt, M) #Unused
    #println("Update")#BUG Debug line



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


    best
end




end
