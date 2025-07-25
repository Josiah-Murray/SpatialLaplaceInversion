using Plots
using BenchmarkTools
Plots.plotlyjs()
M=100

setprecision(Int(ceil(log(2,10)*2.1*M))) #Note: this is not 'thread safe'
tVal = BigFloat(9.9)



#gwr_algorithm = gwr_package
#gwr_algorithm = my_gwr
#gwr_algorithm = gwr_package_arbnumerics
#gwr_algorithm = gwr_package_GThread
#gwr_algorithm = gwr_package_GρThread
gwr_algorithm = gwr_package_iter


#results = @time "Cohen test for $gwr_algorithm" CohenSuite( (LaplaceFunction, t) -> gwr_algorithm(LaplaceFunction,t,M), tVal)

results = @btime CohenSuite( (LaplaceFunction, t) -> gwr_algorithm(LaplaceFunction,t,M), tVal)

errors = abs.(results[:,1] - results[:,2])

p = scatter(errors)
display(p)

for i in 1:35
  println(i,": ", errors[i])
end


##

tNum = 100
tVals = LinRange(0.01, tVal, tNum)
timeTest = [gwr_algorithm(Cohen12, t, M) for t ∈ tVals]
timeTest_exact = [Cohen12_exact(t) for t ∈ tVals]

q = plot(tVals, timeTest, linewidth = 4, box = :on, labels = "Approx.")
q = plot!(tVals, timeTest_exact, linewidth = 3, linestyle = :dash, labels = "Exact")
display(q)
