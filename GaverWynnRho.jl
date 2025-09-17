module GaverWynnRho


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


end
