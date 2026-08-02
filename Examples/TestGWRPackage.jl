using InverseLaplace
using Plots
include("../NILaplace.jl")
using .NILaplace


M = 120
setprecision(3*M)
tVal = BigFloat(9.9)

results = NILaplace.TestFunctions.CohenSuite( (LaplaceFunction, t) -> NILaplace.GaverWynnRho.GWR(LaplaceFunction, t, M; shift_parameter = 1), tVal)

errors = abs.(results[:,1] - results[:,2])./(abs.(results[:,1]))

p = scatter(errors, legend = false)
display(p)

for i in 1:35
  println(i,": ", errors[i])
end



tNum = 100
tVals = LinRange(0.01, tVal, tNum)
ft = t-> NILaplace.GaverWynnRho.GWR(NILaplace.TestFunctions.Cohen12, t, M; shift_parameter = 0.1)
timeTest = real([ft(t) for t ∈ tVals])
timeTest_exact = real([NILaplace.TestFunctions.Cohen12_exact(t) for t ∈ tVals])

q = plot(tVals, timeTest, linewidth = 4, box = :on, labels = "Approx.")
q = plot!(tVals, timeTest_exact, linewidth = 3, linestyle = :dash, labels = "Exact")
display(q)
