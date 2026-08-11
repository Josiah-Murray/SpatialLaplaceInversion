module SpatialLaplaceInversion

using SpecialFunctions


export GWR, GWR_array, DQ, DQ_Array, GenerateWeeksApproximation, EvalWeeks, CohenSuite, CohenSuiteContDiff

include("GaverWynnRho.jl")
include("TestFunctions.jl")
include("Weeks.jl")
include("DirectQuadrature.jl")

using .GaverWynnRho: GWR #Can now access as NILaplace.GWR
using .GaverWynnRho: GWR_array

using .DirectQuadrature: DQ
using .DirectQuadrature: DQ_Array

using .Weeks: GenerateWeeksApproximation
using .Weeks: EvalWeeks

using .TestFunctions: CohenSuite
using .TestFunctions: CohenSuiteContDiff



end
