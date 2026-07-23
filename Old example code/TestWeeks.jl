using Revise
using Plots
include(joinpath(@__DIR__,"../NILaplace.jl"))
using .NILaplace
using InverseLaplace: Weeks

tVals = LinRange(0,1,100)

testFunctionArray = s-> [[NILaplace.TestFunctions.Cohen1(s),
                          NILaplace.TestFunctions.Cohen6(s)];;
                          [NILaplace.TestFunctions.Cohen15(s),
                          NILaplace.TestFunctions.Cohen25(s)]]


exactFunctionArray = t-> [[NILaplace.TestFunctions.Cohen1_exact(t),
                          NILaplace.TestFunctions.Cohen6_exact(t)];;
                          [NILaplace.TestFunctions.Cohen15_exact(t),
                          NILaplace.TestFunctions.Cohen25_exact(t)]]

ApproxStruct = NILaplace.Weeks.GenerateWeeksApproximation(
  testFunctionArray,
  181,
  1/2,
  1
)

#println("Exact: ", NILaplace.TestFunctions.Cohen25_exact(t))
#println("Approx: ", NILaplace.EvalWeeks(ApproxStruct, t))
println("Exact t=0: ", exactFunctionArray(3.2))
println("Approx t=0: ", NILaplace.EvalWeeks(ApproxStruct, 3.2))
