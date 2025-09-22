using .NILaplace: GWR_array
using .NILaplace: GWR

TestFuncArray = s -> [[NILaplace.TestFunctions.Cohen1(s), NILaplace.TestFunctions.Cohen2(s)];; [NILaplace.TestFunctions.Cohen3(s), NILaplace.TestFunctions.Cohen4(s)]]




multDimAnswer = GWR_array(TestFuncArray, 1.0, 20)
println("Multidimensional answer: ", multDimAnswer)
GWR(NILaplace.TestFunctions.Cohen1, 1.0, 20)
