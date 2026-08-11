using Revise
using Plots

using SpatialLaplaceInversion

tVals = LinRange(0,1,100)

testFunctionArray = s-> [[SpatialLaplaceInversion.TestFunctions.Cohen1(s),
                          SpatialLaplaceInversion.TestFunctions.Cohen6(s)];;
                          [SpatialLaplaceInversion.TestFunctions.Cohen15(s),
                          SpatialLaplaceInversion.TestFunctions.Cohen25(s)]]


exactFunctionArray = t-> [[SpatialLaplaceInversion.TestFunctions.Cohen1_exact(t),
                          SpatialLaplaceInversion.TestFunctions.Cohen6_exact(t)];;
                          [SpatialLaplaceInversion.TestFunctions.Cohen15_exact(t),
                          SpatialLaplaceInversion.TestFunctions.Cohen25_exact(t)]]

ApproxStruct = GenerateWeeksApproximation(
  testFunctionArray,
  181,
  1/2,
  1
)

println("Exact t=0: ", exactFunctionArray(3.2))
println("Approx t=0: ", EvalWeeks(ApproxStruct, 3.2))
