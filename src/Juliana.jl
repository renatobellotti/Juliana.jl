module Juliana
export PatientInfo, getOptimisationConfiguration

include("DistanceCalculations.jl")
include("DataLoading.jl")
include("DvhLosses.jl")
include("IdealDoseDistribution.jl")
include("Losses.jl")
include("fiona_standalone/FionaStandalone.jl")
include("OptimisationHelpers.jl")
include("Plotting.jl")

end
