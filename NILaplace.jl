module NILaplace

include("GaverWynnRho.jl")
include("TestFunctions.jl")

using .GaverWynnRho: GWR #Can now access as NILaplace.GWR
using .GaverWynnRho: GWR_array

end
