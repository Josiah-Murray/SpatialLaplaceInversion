using Plots

maxTime = 5
numTVals = 100

tVals = LinRange(maxTime/numTVals,maxTime,numTVals)
M = 32


F(s) = 1/(s+1)
f(t) = exp(-t)

F(s) = 1/((s-1)^2+1)
f(t) = exp(t)sin(t)

exact = [f(t_val) for t_val ∈ tVals ]

approx = [gwr(F, t_val, M) for t_val ∈ tVals]


dif_error = approx - exact
abs_error = abs.(approx-exact)

p_sol = plot(tVals,exact, linewidth = 4, box = :on)
p_sol = plot!(tVals, approx, linewidth = 4, linestyle = :dash, color=:black, ylims = [-1,1]*maximum(abs.(exact)))
display(p_sol)

p_error = plot(tVals, dif_error, linewidth = 4, color = :red, box=:on, ylims = [-maximum(abs_error)*1.2, maximum(abs_error)*1.2])
display(p_error)
