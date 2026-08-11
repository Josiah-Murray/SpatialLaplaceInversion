using Revise
using SpatialLaplaceInversion

t = 2.0

R = 100
N = 1000

TestFunc = s -> SpatialLaplaceInversion.TestFunctions.Cohen3(s)
exact = SpatialLaplaceInversion.TestFunctions.Cohen3_exact(t)



approx = DQ(TestFunc, t, R, N)
println("Approx: ", approx)
println("Exact: ", exact)
