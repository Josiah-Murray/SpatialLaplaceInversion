using Revise
include(joinpath(@__DIR__,"../NILaplace.jl"))
using .NILaplace: DQ

t = 2.0

R = 100
N = 1000

TestFunc = s -> NILaplace.TestFunctions.Cohen3(s)
exact = NILaplace.TestFunctions.Cohen3_exact(t)



approx = DQ(TestFunc, t, R, N)
println("Approx: ", approx)
println("Exact: ", exact)
