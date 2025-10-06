using Revise
using Plots
include(joinpath(@__DIR__,"../NILaplace.jl"))
using .NILaplace
using InverseLaplace: Weeks

tVals = LinRange(0,1,100)

ApproxStruct = NILaplace.WeeksApproximation(
  NILaplace.TestFunctions.Cohen25,
  181,
  1/2,
  1
)

#println("Exact: ", NILaplace.TestFunctions.Cohen25_exact(t))
#println("Approx: ", NILaplace.EvalWeeks(ApproxStruct, t))
println("Exact t=0: ", NILaplace.TestFunctions.Cohen25_exact(0))
println("Approx t=0: ", NILaplace.EvalWeeks(ApproxStruct, 0))
println("Correct approx t=0: ", sum(ApproxStruct.coefficients))

exact = [NILaplace.TestFunctions.Cohen25_exact(t) for t in tVals]
approx = [NILaplace.EvalWeeks(ApproxStruct, t) for t in tVals]

plot(tVals, [exact approx], label = ["Exact" "Approx"])
