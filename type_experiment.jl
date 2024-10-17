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

using LabelledArrays

p = @SLVector [10.0, 28.0, 8 / 3] (:σ, :ρ, :β)

p = @LArray [10.0, 28.0, 8 / 3] (:σ, :ρ, :β)

ps = @SLVector [10.0, 28.0, 8 / 3] (:σ, :ρ, :β)

h = @LArray [10.0, 28.0, 8 / 3]

using StaticArrays

arr = [2.4,3.14,2.71828]

arr = SizedArray(SA[2.4,3.14,2.71828])

arr_2 = MArray(SA[2.4,3.14,2.71828])

using BenchmarkTools

N = 1e3

@btime sum((x -> x^2).(1:N))
@btime mapreduce(x -> x^2,+,1:N)

@code_warntype sum((x -> x^2).(1:N))

@code_llvm sum((x -> x^2).(1:N))
@code_llvm mapreduce(x -> x^2,+,1:N)




B = u0.B
dB = similar(B) / 1u"s"
A = u0.A
function  test_mapreduce(B,d_x,n)
    mapreduce(x -> d_x[x]*(x-1)*B[x],+,1:n+1)
end

function test_sum(B,d_x,n)
    sum((x -> model_params.d_x[x]*(x-1)*B[x]).(1:model_params.n+1))
end

@btime test_mapreduce(B,model_params.d_x,model_params.n)
@btime test_sum(B,model_params.d_x,model_params.n)

function dB_test(B,dB,A,p)
    carrying_coefficient = (p.C - sum(B)) / p.C
    rho_fun = x -> 2*sum(p.f[x:p.n,x] .* B[x:p.n] .* p.r_x[x:p.n]) * carrying_coefficient
     dB_fun = x -> binding_coefficient * (n - (x-1) + 1) * A * B[x-1] -
     p.k_r * (x-1) * B[x] -
     binding_coefficient * (p.n - (x-1)) * A * B[x] +
    p.k_r * ((x-1)+1) * B[x + 1] + rho_fun(x) -
     p.r_x[x]*B[x]*carrying_coefficient-
     p.d_x[x]*B[x]
    bB[1] =  -binding_coefficient * p.n * A * B[1] +
     p.k_r*((1-1)+1)*B[2] +
      rho_fun(1) -
      p.r_x[1]*B[1]*carrying_coefficient - p.d_x[1]*B[1]
    dB[2:p.n] .= dB_fun.(2:p.n)
    bB[p.n+1] = binding_coefficient*A*B[p.n-1] -
     p.k_r * p.n * B[p.n] +
     rho_fun(p.n) -
     p.r_x[p.n]*B[p.n]*carrying_coefficient - p.d_x[p.n]*B[p.n]
end

dB_test(B,dB,model_params)

B[101] = 1e7cell
B[30] = 1e7cell

function rho_fun(x,p,carrying_coefficient)
     2*sum(p.f[x,x:p.n+1] .* B[x:p.n+1] .* p.r_x[x:p.n+1]) *
     carrying_coefficient
end
carrying_coefficient = (model_params.C - sum(B)) / model_params.C
rho_fun(20,model_params,carrying_coefficient)
findlast(model_params.r_x .> zero(model_params.r_x))
model_params.r_x[30] * B[30] * model_params.f[20,30]

B[30]
model_params.r_x[30]
using LinearAlgebra
@btime model_params.f * B ⋅ model_params.r_x

model_params.f[20,30]

model_params.r_x .* B
model_params.f[10:model_params.n+1,10]
model_params.f[begin:end,10]

test_fun_begin = (x,B) -> B[x]

using SpecialFunctions
N=Integer(1e2)
e_mapreduce(x,N) = mapfoldr(n -> x^n/factorial(big(n)),+,0:N)
e_broadcast(x,N) = sum(x.^(0:N)./ factorial.(big.(0:N)))
@btime e_mapreduce(10,N)
@btime e_broadcast(10,N)
@btime exp(10)

include("core_model.jl")
d_u0 = similar(u0) ./ 1u"s"

ode_system!(d_u0, u0, model_params, 0u"s")

@code_warntype ode_system!(d_u0, u0, model_params, 0u"s")
using Cthulhu
@descend ode_system!(d_u0, u0, model_params, 0u"s")