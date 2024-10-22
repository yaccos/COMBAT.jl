using Unitful
import RecursiveArrayTools

# Copy-catted from DiffEqBase DiffEqBaseUnitfulExt.jl
Value(x::Number) = x
Value(x::Type{T}) where {T <:Number} = T
Value(x::Type{Unitful.AbstractQuantity{T, D, U}}) where {T, D, U} = T
Value(x::Unitful.AbstractQuantity) = x.val

# If we do not specify this as a subtype of AbstractVector, the Broadcast machinery will try to convert it into
# a broadcastable representation, but we make this type to be broadcastable as-is due to our customizations of
# copy and copyto!
mutable struct DiscreteSimulationVariables{T<:Number, U<:Number, V<:Number, W<:Number,Z<:Number} <: AbstractVector{T}
    A::U
    T::V
    AT::W
    B::Vector{Z}
    function DiscreteSimulationVariables(args...)
        arg_types = map(RecursiveArrayTools.recursive_bottom_eltype,args)
        T = promote_type(arg_types...)
        new{T,arg_types...}(args...)
    end
end

n_scalars(::DiscreteSimulationVariables{T, U, V, W, Z}) where {T, U, V, W, Z} =  length(fieldnames(DiscreteSimulationVariables{T, U, V, W, Z})) - 1

# Define required AbstractVector methods
Base.length(dsv::DiscreteSimulationVariables) = n_scalars(dsv) + length(dsv.B)
Base.size(dsv::DiscreteSimulationVariables) = (length(dsv),)

Base.firstindex(dsv::DiscreteSimulationVariables) = 1
Base.lastindex(dsv::DiscreteSimulationVariables) = length(dsv)

function Base.getindex(dsv::DiscreteSimulationVariables{T}, idx) where {T}
    if idx == 1
        dsv.A
    elseif idx == 2
        dsv.T
    elseif idx == 3
        dsv.AT
    else
        dsv.B[idx+firstindex(dsv.B)-1-n_scalars(dsv)]
    end
end

@inline Base.@constprop :aggressive function Base.copy(dsv::DiscreteSimulationVariables)
    results = map(sym -> copy(getfield(dsv,sym)),(:A,:T,:AT,:B))
    DiscreteSimulationVariables(results...)
end

@inline Base.@constprop :aggressive function Base.copy!(dst::DiscreteSimulationVariables, src::DiscreteSimulationVariables)
    immutable_syms = (:A,:T,:AT)
    map(immutable_syms) do sym
        res = getfield(src, sym)
        setfield!(dst, sym, res)
    end
    recursivecopy!(dst.B,src.B)
end

Base.setindex!(dsv::DiscreteSimulationVariables{T, U, V, W, Z}, val, idx) where {T, U, V, W, Z} = begin
    if idx == 1
        dsv.A = U(val)
    elseif idx == 2
        dsv.T = V(val)
    elseif idx == 3
        dsv.AT = W(val)
    else
        dsv.B[idx-n_scalars(dsv)] = Z(val)
    end
end


Base.@constprop :aggressive function recursivecopy!(A::DiscreteSimulationVariables,B::DiscreteSimulationVariables)
    immutable_syms = (:A,:T,:AT)
    map(immutable_syms) do sym
        res = getfield(B, sym)
        setfield!(A, sym, res)
    end
    recursivecopy!(A.B,B.B)
end

Base.@constprop :aggressive function recursivecopy(A::DiscreteSimulationVariables) 
    results = map(sym -> recursivecopy(getfield(A,sym)),(:A,:T,:AT,:B))
    DiscreteSimulationVariables(results...)
end



# Implement broadcasting support
# Base.axes(x::DiscreteSimulationVariables) = (Base.OneTo(n_scalars(x)+length(x.B)),)
Base.BroadcastStyle(::Type{<:DiscreteSimulationVariables}) = Broadcast.Style{DiscreteSimulationVariables}()
Base.BroadcastStyle(::Type{<:DiscreteSimulationVariables{T, U, V, W, Z}}) where {T, U, V, W, Z} = Broadcast.Style{DiscreteSimulationVariables{T, U, V, W, Z}}()
Base.BroadcastStyle(::Type{<:DiscreteSimulationVariables{T, U, V, W, Z}}, ::Any) where {T, U, V, W, Z} = Broadcast.Style{DiscreteSimulationVariables{T, U, V, W, Z}}()
Base.BroadcastStyle(::Broadcast.Style{DiscreteSimulationVariables},::Base.Broadcast.BroadcastStyle) = Broadcast.Style{DiscreteSimulationVariables}()
# Base.BroadcastStyle(::Type{<:DiscreteSimulationVariables}) = Base.Broadcast.DefaultArrayStyle{1}()


function Base.similar(x::DiscreteSimulationVariables{T, U, V, W, Z}) where {T, U, V, W, Z}
    DiscreteSimulationVariables(zero(U),zero(V),zero(W),similar(Array{Z}, axes(x.B)))
end

function Base.zero(x::DiscreteSimulationVariables{T, U, V, W, Z}) where {T, U, V, W, Z}
    zero_structures = map(sym -> zero(getfield(x,sym)),(:A,:T,:AT,:B))
    DiscreteSimulationVariables(zero_structures...)
end

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables}})
    # Scan the inputs for the DiscreteSimulationVariables:
    x = find_simulation_variables(bc)
    U = typeof(x.A)
    V = typeof(x.T)
    W = typeof(x.AT)
    Z = eltype(x.B)
    
    DiscreteSimulationVariables(zero(U),zero(V),zero(W),similar(Array{Z}, axes(x.B)))
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
@inline unpack_args(x::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables}}, field::Symbol) = Broadcast.broadcasted(x.f,unpack_args(x.args,field)...)
@inline unpack_args(x::Tuple, field::Symbol) = map(element -> unpack_args(element, field), x)

# The constant propagation is appearently important to ensure type stability
# Otherwise the field symbol does not get propagated, and hence Julia is unable to infer the type returned by getfield

# Using the low-level functions Broadcast.broadcasted or Broadcast.Broadcasted incur considerable
# overhead due to some oddities in the Julia compiler when the arg tuple is not a bitset
# This one is for copies and immutable fields
Base.@constprop :aggressive get_field_res(f, args::Tuple, field::Symbol) = broadcast(f,unpack_args(args, field)...)
# This one is for copying to a mutable field
Base.@constprop :aggressive copyto_field_res!(dest, f, args::Tuple, field::Symbol) = broadcast!(f, dest, unpack_args(args, field)...)



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

@inline Base.@constprop :aggressive function Base.copyto!(dest::DiscreteSimulationVariables,bc::Broadcast.Broadcasted{Broadcast.Style{DiscreteSimulationVariables}, Axes, F, Args}) where {Axes,F,Args<:Tuple}
    f = bc.f
    args = bc.args
    immutable_syms = (:A,:T,:AT)
    # Iterating over immutable symbols which cannot be copied to
    map(immutable_syms) do sym
        res = get_field_res(f,args,sym)
        setfield!(dest, sym, res)
    end
    # Using copto_field_res incurs an overhead, so this is a fine-tuned (albeit not that generalizable) solution for our case
    dest.B .= f.(unpack_args(args,:B)...)
    # copyto_field_res!(dest.B,f, args, :B)
    dest
end

# Adopted from https://scientificcoder.com/user-defined-show-method-in-julia

function Base.show(io::IO, obj::DiscreteSimulationVariables)
    print_object(io, obj, multiline = false)
end

# the 3-argument show used by display(obj) on the REPL
function Base.show(io::IO, mime::MIME"text/plain", obj::DiscreteSimulationVariables)
    # you can add IO options if you want
    multiline = get(io, :multiline, true)
    print_object(io, obj, multiline = multiline)
end

function print_object(io::IO, obj::DiscreteSimulationVariables; multiline::Bool)
    if multiline
        print(io, "Drug effect system with") # or call summary(io, obj)
        print(io, "\n  ")
        print(io, "Antibiotic: $(obj.A)")
        print(io, "\n  ")
        print(io, "Extracellular free target: $(obj.T)")
        print(io, "\n  ")
        print(io, "Extracellular bound target: $(obj.AT)")
        print(io, "\n  ")
        print(io, "Number of bacteria in each compartment:")
        n = 0
        for i in eachindex(obj.B)
            print(io, "\n  ")
            print(io, "B_$(n): $(obj.B[i])")
            n += 1
        end
    else
        # write something short, or go back to default mode
        Base.show_default(io, obj)
    end
end
