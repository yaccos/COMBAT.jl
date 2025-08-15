using Unitful
include("heterogeneous_vector.jl")

x_1 = HeterogeneousVector(u=3.1u"m", v=5.2u"s")
x_2 = HeterogeneousVector(u=1.8u"m", v=8.44u"s")
a = 4.5
b = 9.1
y = zero(x_1)
y .= a .* x_1 .+ b .* x_2
print(y)
