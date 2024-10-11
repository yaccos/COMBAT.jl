using BenchmarkTools
using Unitful

# Base.@constprop :aggressive f(x::AbstractArray,y::AbstractArray) = sum(x) * sum(y)

Base.@constprop :aggressive f(x::AbstractArray,y::AbstractArray) = x[1] * y[1]

Base.@constprop :aggressive f(x::Number,y::Number) = x * y


x_unitful = [1u"m"]

y_unitful =  [3u"m"]


x_unitless =  [1]

y_unitless = [3]


@btime $x_unitful .+ $y_unitful
@btime $x_unitless .+ $y_unitless
