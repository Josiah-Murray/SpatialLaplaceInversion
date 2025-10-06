module NILaplace

include("GaverWynnRho.jl")
include("TestFunctions.jl")
include("Weeks.jl")

using .GaverWynnRho: GWR #Can now access as NILaplace.GWR
using .GaverWynnRho: GWR_array


using .Weeks: WeeksApproximation
using .Weeks: EvalWeeks


end
