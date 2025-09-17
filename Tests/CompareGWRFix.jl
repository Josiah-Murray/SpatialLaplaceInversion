using Plots
using BenchmarkTools
Plots.plotlyjs()

MVals = 10:2:64

Better = Array{Int}(undef, length(MVals))
Same = Array{Int}(undef, length(MVals))
Worse = Array{Int}(undef, length(MVals))



for Mi in eachindex(MVals)
  M = MVals[Mi]

  setprecision(Int(ceil(log(2,10)*2.1*M))) #Note: this is not 'thread safe'
  tVal = BigFloat(15)



  gwr_algorithm = gwr_package_GFFix

  results = CohenSuite( (LaplaceFunction, t) -> gwr_algorithm(LaplaceFunction,t,M), tVal)

  errors = abs.(results[:,1] - results[:,2])



  #OLD

  results_old = CohenSuite( (LaplaceFunction, t) -> gwr_package(LaplaceFunction,t,M), tVal)

  errors_old = abs.(results_old[:,1] - results_old[:,2])

  ##

  count_better = 0
  count_same = 0
  count_worse = 0
  for i in eachindex(errors)
    if errors[i] > errors_old[i]
       count_worse += 1
    elseif errors[i] < errors_old[i]
       count_better += 1
    else
       count_same += 1
    end
  end

  Better[Mi] = count_better
  Same[Mi] = count_same
  Worse[Mi] = count_worse


end


println("---------------------------------------")
println("   M: |  Better:  |  Same:  |  Worse:  ")
println("------|-----------|---------|----------")
for i in eachindex(MVals)
  println( "  ",  length(string(MVals[i]))== 1 ? " " : "", MVals[i],"  |    ",
    length(string(Better[i]))== 1 ? " " : "", Better[i], "     |   " ,
    length(string(Same[i]))== 1 ? " " : "", Same[i], "    |     ",
    length(string(Worse[i]))== 1 ? " " : "", Worse[i])
end
