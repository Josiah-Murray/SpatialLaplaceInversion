using Revise
using SpatialLaplaceInversion
t = 2.0


TestFuncArray = s -> [[SpatialLaplaceInversion.TestFunctions.Cohen1(s), SpatialLaplaceInversion.TestFunctions.Cohen2(s)];; [SpatialLaplaceInversion.TestFunctions.Cohen3(s), SpatialLaplaceInversion.TestFunctions.Cohen4(s)]]
exact = [[SpatialLaplaceInversion.TestFunctions.Cohen1_exact(t), SpatialLaplaceInversion.TestFunctions.Cohen2_exact(t)];;  [SpatialLaplaceInversion.TestFunctions.Cohen3_exact(t), SpatialLaplaceInversion.TestFunctions.Cohen4_exact(t) ]]

R = 100
N = 20000


multDimAnswer = real.(DQ_Array(TestFuncArray, t, R, N, γ = 0.01))
println("Multidimensional answer: ", multDimAnswer)
println("Exact: ", exact)
