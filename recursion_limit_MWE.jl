using BenchmarkTools

struct MyType{T} <: AbstractVector{Number}
    val::T
end

# Only intended for displaying the type
Base.getindex(x::MyType, idx) = x.val[idx]

Base.BroadcastStyle(::Type{<:MyType}) = Broadcast.Style{MyType}()


@inline Base.length(x::MyType) = length(x.val)
Base.size(x::MyType) = (length(x),)

@inline unpack_args(x::MyType) = x.val
@inline unpack_args(x::Broadcast.Broadcasted{Broadcast.Style{MyType}}) = Broadcast.broadcasted(x.f,unpack_args(x.args)...)
@inline unpack_args(x::Tuple) = map(element -> unpack_args(element), x)


@inline function Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{MyType}, Axes, F, Args}) where {Axes,F,Args<:Tuple}
    f = bc.f
    args = bc.args
    new_val = f.(unpack_args(args)...)
    MyType(new_val)
end

@inline function Base.copyto!(dest::MyType,bc::Broadcast.Broadcasted{Broadcast.Style{MyType}, Axes, F, Args}) where {Axes,F,Args<:Tuple}
    f = bc.f
    args = bc.args
    dest.val .= f.(unpack_args(args)...)
    dest
end

x = MyType([1.0, 2.0, 3.0])

function f(val, N)
    for i in 1:N
        val .+= val
    end
end

N = 1e9

@btime f(x, N)
@btime f(x.val, N)

@code_warntype f(x, N)

@profview f(x, N)
@profview f(x.val, N)

using Cthulhu
@descend f(x, N)
@descend f(x.val, N)
@descend broadcast!(+,x,x,x)

broadcast!(+,x,x,x)
