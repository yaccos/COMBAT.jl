using Unitful
using MultiScaleArrays

mutable struct DiscreteSimulationVariables{T<:Number, U<:Number, V<:Number, W<:Number} <: AbstractMultiScaleArrayLeaf{Float64}
    A::T
    T::U
    AT::V
    B::Vector{W}
end

function create_simulation_variables()
    A = 1.0u"m"   # meters
    T = 2.0u"s"   # seconds
    AT = 3.0u"m/s"  # meters per second
    B = [1.0u"m", 2.0u"m", 3.0u"m"]  # vector of meters

    return DiscreteSimulationVariables(A, T, AT, B)
end

n_scalars = length(fieldnames(DiscreteSimulationVariables)) - 1

# Define required AbstractVector methods
Base.length(x::DiscreteSimulationVariables) = n_scalars + length(x.B)
Base.size(x::DiscreteSimulationVariables) = size(length(x),)
function Base.similar(x::DiscreteSimulationVariables{T, U, V, W}) where {T, U, V, W}
    DiscreteSimulationVariables(T(0),U(0),V(0),similar(Array{W}, axes(x.B)))
end

function Base.similar(x::DiscreteSimulationVariables{T, U, V, W}, ::Type{ElType}) where {T, U, V, W, ElType}
    zero_value = zero(ElType)
    DiscreteSimulationVariables(T(zero_value),U(zero_value),V(zero_value),similar(Array{W}, axes(x.B)))
end

Base.getindex(dsv::DiscreteSimulationVariables, idx) = dsv.B[idx]
Base.setindex!(dsv::DiscreteSimulationVariables, val, idx) = (dsv.B[idx] = val)

# Example usage with specific Unitful types

measurements = [1.0u"m",2.0u"s",3.0u"m/s", 1.0u"m", 2.0u"m", 3.0u"m"]

import Random
rng = Random.MersenneTwister(3052)
n_numbers = Int64(1e6)
random_numbers = exp.(randn(rng,Float64,(n_numbers,)))
f(x,random_numbers) = begin
    res = deepcopy(x)
    for i in random_numbers
        res .*= i
    end
    res
end

@time f(measurements,random_numbers)

measurements_raw = ustrip.(measurements)

@time f(measurements_raw,random_numbers)

measurements_homogenuous = measurements_raw * 1.0u"m"

@time f(measurements_homogenuous,random_numbers)

l=Broadcast.Broadcasted(+,(val,5))

fieldnames(l)

l.f

l.args

l.axes

fieldnames(l)


Broadcast.Broadcasted

