using Revise
include(joinpath(@__DIR__,"../NILaplace.jl"))
using .NILaplace
t = 2.0


TestFuncArray = s -> [[NILaplace.TestFunctions.Cohen1(s), NILaplace.TestFunctions.Cohen2(s)];; [NILaplace.TestFunctions.Cohen3(s), NILaplace.TestFunctions.Cohen4(s)]]
exact = [[NILaplace.TestFunctions.Cohen1_exact(t), NILaplace.TestFunctions.Cohen2_exact(t)];;  [NILaplace.TestFunctions.Cohen3_exact(t), NILaplace.TestFunctions.Cohen4_exact(t) ]]

R = 100
N = 20000


multDimAnswer = real.(NILaplace.DQ_Array(TestFuncArray, t, R, N, γ = 0.01))
println("Multidimensional answer: ", multDimAnswer)
println("Exact: ", exact)
