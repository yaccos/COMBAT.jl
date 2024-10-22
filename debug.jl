include("array_interface.jl")
import RecursiveArrayTools
function create_simulation_variables()
    A = 1.0u"m"   # meters
    T = 2.0u"s"   # seconds
    AT = 3.0u"m/s"  # meters per second
    B = [1.0u"m", 2.0u"m", 3.0u"m"]  # vector of meters

    return DiscreteSimulationVariables(A, T, AT, B)
end

val = create_simulation_variables()
val_copy = RecursiveArrayTools.recursivecopy(val)
