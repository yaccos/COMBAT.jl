using Unitful
using MultiScaleArrays


function create_simulation_variables()
    A = 1.0u"m"   # meters
    T = 2.0u"s"   # seconds
    AT = 3.0u"m/s"  # meters per second
    B = [1.0u"m", 2.0u"m", 3.0u"m"]  # vector of meters

    return DiscreteSimulationVariables(A, T, AT, B)
end

val = create_simulation_variables()
val .+= val



@code_warntype val .+ val

@code_warntype val .* val

val_2 = copy(val)

val_3 = copy(val)

@code_warntype copy(val)

@code_warntype broadcast!(*,val_3, val_2, val_2)


using Cthulhu

Cthulhu.@descend broadcast(+, val,val)

Cthulhu.@descend broadcast!(+,val,val,val)

Cthulhu.@descend broadcast!(+,val_3, val_2, val_2)

@time val_3 .= val_2 .+ val_2

@time val_4 = val_2 .+ val_2

@code_warntype broadcast!(+,val_3, val_2, val_2)

@code_warntype fieldnames(val)

#Cthulhu.@descend fieldnames(val)

#Cthulhu.@descend n_scalars(val)

typeof(val)

A = rand(100)

B = zeros(100)

C = copyto!(B,A)

map(i -> setindex!(A,i,i), eachindex(A))

A

typeof(typeof(val))

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

@time f(val,N)

@time g(val,N)

@profview_allocs g(val,N)

@profview_allocs f(val,N)

@profview g(val,N)

@profview g(val,N)

f(x) = get_field_res(+,(x,x),:A)

h(x) = get_field_res(+,(x,x),:B)

function prof_fun(x,sym, N)
    for i in 1:N
        get_field_res(+,(x,x),sym)
    end
end

@time prof_fun(val, :B, N)

@time f(val)
@time h(val)
@time get_field_res(+,(val,val),:A)

map(x -> x^2, (5,9,19))

map(x -> x^2, ())
