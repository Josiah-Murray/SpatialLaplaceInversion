using Revise
include(joinpath(@__DIR__,"../NILaplace.jl"))
using .NILaplace: GWR_array
using .NILaplace: GWR

t = 2.0


TestFuncArray = s -> [[NILaplace.TestFunctions.Cohen1(s), NILaplace.TestFunctions.Cohen2(s)];; [NILaplace.TestFunctions.Cohen3(s), NILaplace.TestFunctions.Cohen4(s)]]
exact = [[NILaplace.TestFunctions.Cohen1_exact(t), NILaplace.TestFunctions.Cohen2_exact(t)];;  [NILaplace.TestFunctions.Cohen3_exact(t), NILaplace.TestFunctions.Cohen4_exact(t) ]]



multDimAnswer = GWR_array(TestFuncArray, t, 20)
println("Multidimensional answer: ", multDimAnswer)
GWR(NILaplace.TestFunctions.Cohen1, 1.0, 20)
