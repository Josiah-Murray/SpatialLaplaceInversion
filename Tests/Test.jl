using Plots
M=40

tVal = 9.9

results = CohenSuiteContDiff( (LaplaceFunction, t) -> my_gwr(LaplaceFunction,t,M), tVal)

errors = abs.(results[:,1] - results[:,2])

p = bar(eachindex(errors), errors)
p = annotate!(eachindex(errors), string.(results[:,3]), :bottom)
display(p)

for i in eachindex(errors)
  println(results[i,3],": ", errors[i])
end


##
tNum = 100
tVals = LinRange(0.01, tVal, tNum)
timeTest = [my_gwr(Cohen19, t, M) for t ∈ tVals]
timeTest_exact = [Cohen19_exact(t) for t ∈ tVals]

q = plot(tVals, timeTest, linewidth = 4, box = :on, labels = "Approx.")
q = plot!(tVals, timeTest_exact, linewidth = 3, linestyle = :dash, labels = "Exact")
display(q)
