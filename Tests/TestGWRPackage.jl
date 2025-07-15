using InverseLaplace
using Plots

M=100

tVal = 9.9

results = CohenSuite( (LaplaceFunction, t) -> GWR(LaplaceFunction,M)(t), tVal)

errors = abs.(results[:,1] - results[:,2])

p = scatter(errors)
display(p)

for i in 1:35
  println(i,": ", errors[i])
end



tNum = 100
tVals = LinRange(0.01, tVal, tNum)
ft = GWR(Cohen12, M)
timeTest = [ft(t) for t ∈ tVals]
timeTest_exact = [Cohen12_exact(t) for t ∈ tVals]

q = plot(tVals, timeTest, linewidth = 4, box = :on, labels = "Approx.")
q = plot!(tVals, timeTest_exact, linewidth = 3, linestyle = :dash, labels = "Exact")
display(q)