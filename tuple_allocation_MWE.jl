using StaticArrays
using BenchmarkTools
using Unitful

struct MyType{T <: AbstractVector}
    val::T
end


f(x::AbstractArray,y::AbstractArray) = sum(x) * sum(y)

@noinline f_wrapper(args) = f(args...)

function f_wrapper_1(x::MyType,y::MyType)
    f(x.val,y.val)
end

function f_wrapper_2(x::MyType,y::MyType)
    args = (x.val,y.val)
    f_wrapper(args)
end

x_static = MyType(@SVector [1u"m",2u"m"])

y_static = MyType(@SVector [3u"m",4u"m"])

x_nonstatic = MyType([1u"m",2u"m"])

y_nonstatic = MyType([3u"m",4u"m"])

@btime f_wrapper_1(x_static,y_static)

@btime f_wrapper_1(x_nonstatic,y_nonstatic)

@btime f_wrapper_2(x_static,y_static)

@btime f_wrapper_2(x_nonstatic,y_nonstatic)

@code_llvm f_wrapper_1(x_nonstatic,y_nonstatic)

@code_llvm f_wrapper_2(x_nonstatic,y_nonstatic)

@code_warntype f_wrapper_1(x_static,y_static) 
