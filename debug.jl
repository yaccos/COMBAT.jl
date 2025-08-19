using Unitful
using BenchmarkTools
include("heterogeneous_vector.jl")

x_1 = HeterogeneousVector(u=3.1u"m", v=5.2u"s")
x_2 = HeterogeneousVector(u=1.8u"m", v=8.44u"s")
a = 4.5
b = 9.1
y_2 = a .* x_1 .+ b .* x_2
y_1 = zero(x_1)
y_1 .= a .* x_1 .+ b .* x_2

function f(val, N)
    for i in 1:N
        val .+= val
    end
end

N = 1e9

f(y_1, N)

# @profview f(y_1, N) 
