module SpatialLaplaceInversion

using SpecialFunctions


export GWR, GWR_array, DQ, DQ_Array, WeeksApproximation, EvalWeeks

include("GaverWynnRho.jl")
include("TestFunctions.jl")
include("Weeks.jl")
include("DirectQuadrature.jl")

using .GaverWynnRho: GWR #Can now access as NILaplace.GWR
using .GaverWynnRho: GWR_array

using .DirectQuadrature: DQ
using .DirectQuadrature: DQ_Array

using .Weeks: WeeksApproximation
using .Weeks: EvalWeeks

#TODO: Export test functions.




end
