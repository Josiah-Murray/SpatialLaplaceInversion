using Plots
<<<<<<< HEAD
M=40
=======
using BenchmarkTools
Plots.plotlyjs()
>>>>>>> GWR-experiments

M = 14

<<<<<<< HEAD
results = CohenSuiteContDiff( (LaplaceFunction, t) -> my_gwr(LaplaceFunction,t,M), tVal)
=======
setprecision(Int(ceil(log(2,10)*2.1*M))) #Note: this is not 'thread safe'
tVal = BigFloat(9.9)



#wr_algorithm = gwr_package
#gwr_algorithm = my_gwr
#gwr_algorithm = gwr_package_arbnumerics
#gwr_algorithm = gwr_package_GThread
#gwr_algorithm = gwr_package_GρThread
#gwr_algorithm = gwr_package_iter
#gwr_algorithm = gwr_package_iter_threads
gwr_algorithm = gwr_package_GFFix

#results = @time "Cohen test for $gwr_algorithm" CohenSuite( (LaplaceFunction, t) -> gwr_algorithm(LaplaceFunction,t,M), tVal)

results = CohenSuite( (LaplaceFunction, t) -> gwr_algorithm(LaplaceFunction,t,M), tVal)
>>>>>>> GWR-experiments

errors = abs.(results[:,1] - results[:,2])

p = bar(eachindex(errors), errors)
p = annotate!(eachindex(errors), string.(results[:,3]), :bottom)
display(p)

<<<<<<< HEAD
for i in eachindex(errors)
  println(results[i,3],": ", errors[i])
end
=======
>>>>>>> GWR-experiments

#OLD

results_old = CohenSuite( (LaplaceFunction, t) -> gwr_package(LaplaceFunction,t,M), tVal)

errors_old = abs.(results_old[:,1] - results_old[:,2])

temp = abs.(errors_old)-abs.(errors)
temp = (temp)./(abs.(errors_old)+abs.(errors))
temp_plot = scatter(temp)
plot!(title ="Error comparison for M = $M",
  legend = false,
  box = :origin,
  ylims = [-1.1,1.1],
  xticks = 1:35
)
display(temp_plot)

##

for i in 1:35
  println(i,": ", errors[i], "  |  ", errors_old[i])
end

count_better = 0
count_same = 0
count_worse = 0
for i in eachindex(errors)
  if errors[i] > errors_old[i]
    global count_worse += 1
  elseif errors[i] < errors_old[i]
    global count_better += 1
  else
    global count_same += 1
  end
end

println("Number better: ", count_better)
println("Number worse: ", count_worse)
println("Number same: ", count_same)



#=
tNum = 100
tVals = LinRange(0.01, tVal, tNum)
<<<<<<< HEAD
timeTest = [my_gwr(Cohen19, t, M) for t ∈ tVals]
timeTest_exact = [Cohen19_exact(t) for t ∈ tVals]
=======
timeTest = [gwr_algorithm(Cohen30, t, M) for t ∈ tVals]
timeTest_exact = [Cohen30_exact(t) for t ∈ tVals]
>>>>>>> GWR-experiments

q = plot(tVals, timeTest, linewidth = 4, box = :on, labels = "Approx.")
q = plot!(tVals, timeTest_exact, linewidth = 3, linestyle = :dash, labels = "Exact")
display(q)
=#
