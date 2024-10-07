using Unitful

mutable struct DiscreteSimulationVariables{T<:Number, U<:Number, V<:Number, W<:Number} <: AbstractVector{Number}
    A::T
    T::U
    AT::V
    B::Vector{W}
end

n_scalars(::DiscreteSimulationVariables) =  length(fieldnames(DiscreteSimulationVariables)) - 1

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
Base.axes(x::DiscreteSimulationVariables) = (Base.OneTo(n_scalars(x)+length(x.B)),)
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

unpack_args(x::Any,::Symbol) = x
unpack_args(::Tuple{}, ::Symbol) = Tuple{}()
# We can assume that we have at least one element
unpack_args(x::Tuple, field::Symbol) = (unpack_args(x[1], field), unpack_args(Base.tail(x), field)...)
unpack_args(x::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables}}, field::Symbol) = get_field_res(x.f,x.args, field)
# This must be a scalar expression in order to be compatible, so there is no harm in materializing right away
unpack_args(x::Broadcast.Broadcasted, ::Symbol) = Base.materialize(x)
unpack_args(x::DiscreteSimulationVariables,field::Symbol) = getfield(x, field)

get_field_res(f, args::Tuple, field::Symbol) = copy(Broadcast.Broadcasted(f,unpack_args(args, field)))



function Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables{T, U, V, W}}, Axes, F, Args}) where {Axes,F,Args<:Tuple,T, U, V, W}
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
    A_res = get_field_res(f,args, :A)
    T_res = get_field_res(f,args, :T)
    AT_res = get_field_res(f,args, :AT)
    B_res = get_field_res(f,args, :B)
    DiscreteSimulationVariables(A_res,T_res,AT_res,B_res)
end

"`x = find_simulation_variables(xs)` returns the first ArrayAndChar among the arguments."


# function Base.similar(bc::Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}}, ::Type{ElType}) where ElType
#     return Vector{ElType}(undef, Base.Broadcast.length(bc))
# end


function create_simulation_variables()
    A = 1.0u"m"   # meters
    T = 2.0u"s"   # seconds
    AT = 3.0u"m/s"  # meters per second
    B = [1.0u"m", 2.0u"m", 3.0u"m"]  # vector of meters

    return DiscreteSimulationVariables(A, T, AT, B)
end
