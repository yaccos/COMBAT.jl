using Unitful

mutable struct DiscreteSimulationVariables{T<:Number, U<:Number, V<:Number, W<:Number} <: AbstractVector{Number}
    A::T
    T::U
    AT::V
    B::Vector{W}
end

n_scalars(::DiscreteSimulationVariables{T, U, V, W}) where {T, U, V, W} =  length(fieldnames(DiscreteSimulationVariables{T, U, V, W})) - 1

# Define required AbstractVector methods
Base.length(dsv::DiscreteSimulationVariables) = n_scalars(dsv) + length(dsv.B)
Base.size(dsv::DiscreteSimulationVariables) = (length(dsv),)
Base.getindex(dsv::DiscreteSimulationVariables, idx) = begin
    if idx == 1
        dsv.A
    elseif idx == 2
        dsv.T
    elseif idx == 3
        dsv.AT
    else
        dsv.B[idx-n_scalars(dsv)]
    end
end

Base.setindex!(dsv::DiscreteSimulationVariables, val, idx) = begin
    if idx == 1
        dsv.A = val
    elseif idx == 2
        dsv.T = val
    elseif idx == 3
        dsv.AT = val
    else
        dsv.B[idx-n_scalars(dsv)] = val
    end
end



# Implement broadcasting support
# Base.axes(x::DiscreteSimulationVariables) = (Base.OneTo(n_scalars(x)+length(x.B)),)
Base.BroadcastStyle(::Type{<:DiscreteSimulationVariables}) = Broadcast.Style{DiscreteSimulationVariables}()
Base.BroadcastStyle(::Type{<:DiscreteSimulationVariables{T, U, V, W}}) where {T, U, V, W} = Broadcast.Style{DiscreteSimulationVariables{T, U, V, W}}()
Base.BroadcastStyle(::Type{<:DiscreteSimulationVariables{T, U, V, W}}, ::Any) where {T, U, V, W} = Broadcast.Style{DiscreteSimulationVariables{T, U, V, W}}()
Broadcast.Style{DiscreteSimulationVariables}()
Base.BroadcastStyle(::Broadcast.Style{DiscreteSimulationVariables},::Base.Broadcast.BroadcastStyle) = Broadcast.Style{DiscreteSimulationVariables}()
# Base.BroadcastStyle(::Type{<:DiscreteSimulationVariables}) = Base.Broadcast.DefaultArrayStyle{1}()


function Base.similar(x::DiscreteSimulationVariables{T, U, V, W}) where {T, U, V, W}
    DiscreteSimulationVariables(T(0),U(0),V(0),similar(Array{W}, axes(x.B)))
end

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables}})
    # Scan the inputs for the DiscreteSimulationVariables:
    x = find_simulation_variables(bc)
    T = typeof(x.A)
    U = typeof(x.T)
    V = typeof(x.AT)
    W = eltype(x.B)
    
    DiscreteSimulationVariables(T(0),U(0),V(0),similar(Array{W}, axes(x.B)))
end



function Base.similar(bc::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables{T, U, V, W}}}) where {T, U, V, W}
    # Scan the inputs for the DiscreteSimulationVariables:
    x = find_simulation_variables(bc)
    DiscreteSimulationVariables(T(0),U(0),V(0),similar(Array{W}, axes(x.B)))
end

# Inspired by MultiScaleArrays.jl
find_simulation_variables(bc::Base.Broadcast.Broadcasted) = find_simulation_variables(bc.args)
find_simulation_variables(args::Tuple) = find_simulation_variables(find_simulation_variables(args[1]), Base.tail(args))
find_simulation_variables(x) = x
find_simulation_variables(::Tuple{}) = nothing
find_simulation_variables(x::DiscreteSimulationVariables, rest) = x
find_simulation_variables(::Any, rest) = find_simulation_variables(rest)

@inline unpack_args(x::Any,::Symbol) = x
@inline unpack_args(x::DiscreteSimulationVariables,field::Symbol) = getfield(x, field)
@inline unpack_args(x::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables}}, field::Symbol) = Broadcast.Broadcasted(x.f,unpack_args(x.args,field))
@inline unpack_args(x::Tuple, field::Symbol) = map(element -> unpack_args(element, field), x)

# The constant propagation is appearently important to ensure type stability
# Otherwise the field symbol does not get propagated, and hence Julia is unable to infer the type returned by getfield 
# This one is for copies and immutable fields
Base.@constprop :aggressive get_field_res(f, args::Tuple, field::Symbol) = copy(Broadcast.Broadcasted(f,unpack_args(args, field)))
# This one is for copying to a mutable field
Base.@constprop :aggressive function copyto_field_res!(dest, f, args::Tuple, field::Symbol)
    clean_args = unpack_args(args,field)
    broadcast = Broadcast.Broadcasted(f,clean_args)
    copyto!(dest, broadcast)
end 



@inline function Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables}, Axes, F, Args}) where {Axes,F,Args<:Tuple}
    f = bc.f
    args = bc.args
    res_args = map(sym -> get_field_res(f,args, sym),(:A,:T,:AT,:B))
    # A_res = get_field_res(f,args, :A)
    # T_res = get_field_res(f,args, :T)
    # AT_res = get_field_res(f,args, :AT)
    # B_res = get_field_res(f,args, :B)
    DiscreteSimulationVariables(res_args...)
end

function Base.copyto!(dest::DiscreteSimulationVariables,bc::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables}, Axes, F, Args}) where {Axes,F,Args<:Tuple}
    f = bc.f
    args = bc.args
    immutable_syms = (:A,:T,:AT)
    # Iterating over immutable symbols which cannot be copied to
    map(immutable_syms) do sym 
        res = get_field_res(f,args,sym)
        setfield!(dest, sym, res)
    end
    copyto_field_res!(dest.B,f, args, :B)
    dest
end



# function Base.similar(bc::Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}}, ::Type{ElType}) where ElType
#     return Vector{ElType}(undef, Base.Broadcast.length(bc))
# end


