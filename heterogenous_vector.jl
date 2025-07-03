using Unitful
using Revise
import RecursiveArrayTools
import Base.zero


# Copy-catted from DiffEqBase DiffEqBaseUnitfulExt.jl
Value(x::Number) = x
Value(x::Type{T}) where {T <:Number} = T
Value(x::Type{Unitful.AbstractQuantity{T, D, U}}) where {T, D, U} = T
Value(x::Unitful.AbstractQuantity) = x.val


# If we do not specify this as a subtype of AbstractVector, the Broadcast machinery will try to convert it into
# a broadcastable representation, but we make this type to be broadcastable as-is due to our customizations of
# copy and copyto!

# Helper function to wrap scalars in Ref for mutability
_make_mutable(x::AbstractArray) = x
_make_mutable(x::Ref) = x
_make_mutable(x) = Ref(x)
_unwrap(x::Ref) = x[]
_unwrap(x) = x
_set_value!(x::Ref, val) = (x[] = val)
_set_value!(x::AbstractArray, val::AbstractArray) = copy!(x,val)
_set_value!(x::AbstractArray, val, idx) = (x[idx] = val)

struct HeterogenousVector{T, S <: NamedTuple} <: AbstractVector{T}
    x::S
    function HeterogenousVector(x::NamedTuple)
        # Wrap scalar fields in Ref for mutability
        mutable_x = map(_make_mutable, x)
        arg_types = map(field -> _unwrap(field) |> RecursiveArrayTools.recursive_bottom_eltype, values(mutable_x))
        T = promote_type(arg_types...)
        new{T, typeof(mutable_x)}(mutable_x)
    end

    # Constructor with keyword arguments
    function HeterogenousVector(; kwargs...)
        x = NamedTuple(kwargs)
        HeterogenousVector(x)
    end

    # Constructor with positional arguments
    function HeterogenousVector(args...)
        names = ntuple(i -> Symbol("field_$i"), length(args))
        x = NamedTuple{names}(args)
        HeterogenousVector(x)
    end
end

function NamedTuple(hv::HeterogenousVector)
    return getfield(hv, :x)
end

# Custom property access for clean external interface
# Note: For accessing the named tuple field x, we must use getfield or invoking NamedTuple
function Base.getproperty(hv::HeterogenousVector, name::Symbol)
    if haskey(NamedTuple(hv), name)
        #  Access field and unwrap if it's a Ref
        field = getfield(NamedTuple(hv), name)
        return _unwrap(field)
    else
        error("HeterogenousVector has no field $name")
    end
end

function Base.setproperty!(hv::HeterogenousVector, name::Symbol, value)
    if haskey(NamedTuple(hv), name)
        # Set field value, wrapping in Ref if it's a scalar
        field = getfield(NamedTuple(hv), name)
        _set_value!(field, value)
    else
        error("HeterogenousVector has no field $name")
    end
end


function Base.propertynames(hv::HeterogenousVector)
    return keys(NamedTuple(hv))
end

Base.pairs(hv::HeterogenousVector{T, S}) where {S, T} = pairs(NamedTuple(hv))

function Base.getindex(hv::HeterogenousVector{T, S}, idx::Int) where {T, S}
    current_idx = 1
    for (name, field) in pairs(hv)
        unwrapped_field = _unwrap(field)
        if unwrapped_field isa AbstractArray
            field_length = length(unwrapped_field)
            if current_idx <= idx < current_idx + field_length
                return unwrapped_field[idx - current_idx + 1]
            end
            current_idx += field_length
        else
            if idx == current_idx
                return unwrapped_field
            end
            current_idx += 1
        end
    end
    throw(BoundsError(hv, idx))
end

function Base.setindex!(hv::HeterogenousVector{T, S}, val, idx::Int) where {T, S}
    current_idx = 1
    for (name, field) in pairs(hv)
        if field isa AbstractArray
            field_length = length(field)
            if current_idx <= idx < current_idx + field_length
                field[idx - current_idx + 1] = val
                return val
            end
            current_idx += field_length
        else
            if idx == current_idx
                # Mutate the Ref directly
                _set_value!(field, val)
                return val
            end
            current_idx += 1
        end
    end
    throw(BoundsError(hv, idx))
end


_field_length(field::Ref) = 1
_field_length(field::AbstractArray) = length(field)
# Update length calculation
Base.length(hv::HeterogenousVector) = sum(_field_length, NamedTuple(hv))

Base.size(hv::HeterogenousVector) = (length(hv),)
Base.firstindex(hv::HeterogenousVector) = 1
Base.lastindex(hv::HeterogenousVector) = length(hv)

_copy_field(field::Ref) = _unwrap(field) |> Ref
_copy_field(field::AbstractArray) = copy(field)
function Base.copy(hv::HeterogenousVector)
    copied_x = map(_copy_field, NamedTuple(hv))
    HeterogenousVector(copied_x)
end



_copy_field!(dst::Ref, src::Ref) = _set_value!(dst, _unwrap(src))
_copy_field!(dst::AbstractArray, src::AbstractArray) = copy!(dst, src)
function Base.copyto!(dst::HeterogenousVector, src::HeterogenousVector)
    # Ensure both have the same structure
    if propertynames(dst) != propertynames(src)
        throw(ArgumentError("HeterogenousVectors must have the same field names"))
    end
    
    for name in propertynames(dst)
        src_field = getfield(NamedTuple(src), name)
        dst_field = getfield(NamedTuple(dst), name)
        
        _copy_field!(dst_field, src_field)
    end
    return dst
end

Base.copy!(dst::HeterogenousVector, src::HeterogenousVector) = Base.copyto!(dst, src)

Base.zero(x::Ref) = _unwrap(x) |> zero
Base.similar(x::Ref) = zero(x)
Base.similar(x::Ref,::Type{ElType}) where {ElType} = Ref(zero(ElType))

function Base.similar(hv::HeterogenousVector{T}) where {T}
    similar_x = map(similar, NamedTuple(hv))
    HeterogenousVector(similar_x)
end

function Base.zero(hv::HeterogenousVector)
    zero_x = map(zero, NamedTuple(hv))
    HeterogenousVector(zero_x)
end

# Broadcasting support for HeterogenousVector
Base.BroadcastStyle(::Type{<:HeterogenousVector{T, S}}) where {T, S} = Broadcast.Style{HeterogenousVector{fieldnames(S)}}()

function Base.BroadcastStyle(::Broadcast.Style{HeterogenousVector{Names1}}, ::Broadcast.Style{HeterogenousVector{Names2}}) where {Names1, Names2} 
    error("Cannot broadcast HeterogenousVectors with different field names: $(Names1) vs $(Names2)")
end

function Base.BroadcastStyle(::Broadcast.Style{HeterogenousVector{Names}}, ::Broadcast.Style{HeterogenousVector{Names}}) where {Names} 
    Broadcast.Style{HeterogenousVector{Names}}()
end

# HeterogenousVector style takes precedence over other broadcast styles
function Base.BroadcastStyle(::Broadcast.Style{HeterogenousVector{Names}}, ::Base.Broadcast.BroadcastStyle) where {Names}
    Broadcast.Style{HeterogenousVector{Names}}()
end

# Helper function to find HeterogenousVector in broadcast arguments
find_heterogenous_vector(bc::Base.Broadcast.Broadcasted) = find_heterogenous_vector(bc.args)
find_heterogenous_vector(args::Tuple) = find_heterogenous_vector(find_heterogenous_vector(args[1]), Base.tail(args))
find_heterogenous_vector(x::Base.Broadcast.Extruded) = x.x
find_heterogenous_vector(x) = x
find_heterogenous_vector(::Tuple{}) = nothing
find_heterogenous_vector(x::HeterogenousVector, rest) = x
find_heterogenous_vector(::Any, rest) = find_heterogenous_vector(rest)

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector{Names}}}) where {Names}
    hv = find_heterogenous_vector(bc)
    similar(hv)
end

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector{Names}}}, ::Type{ElType}) where {Names, ElType}
    hv = find_heterogenous_vector(bc)
    similar_x = map(NamedTuple(hv)) do field
        similar(field, ElType)
    end
    HeterogenousVector(similar_x)
end


# Using generated functions here really helps performance
@generated unpack_args(x::Any,::Symbol) = :(x)
@generated unpack_args(x::HeterogenousVector, field::Symbol) = :(getproperty(x, field))
@generated unpack_args(x::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector{Names}}}, field::Symbol) where {Names} = :(Broadcast.broadcasted(x.f,unpack_args(x.args,field)...))
@generated unpack_args(::Tuple{}, ::Symbol) = :()
@generated unpack_args(x::Tuple, field::Symbol) = :(unpack_args(x[1],field),unpack_args(Base.tail(x),field)...)

# The constant propagation is appearently important to ensure type stability
# Otherwise the field symbol does not get propagated, and hence Julia is unable to infer the type returned by getfield

# Using the low-level functions Broadcast.broadcasted or Broadcast.Broadcasted incur considerable
# overhead due to some oddities in the Julia compiler when the arg tuple is not a bitset
# This one is for copies and immutable fields
@generated get_field_res(f, args::Tuple, field::Symbol) = :(broadcast(f,unpack_args(args, field)...))
# This one is for copying to a mutable field
@generated copyto_field_res!(dest::AbstractArray, f, args::Tuple, field::Symbol) = :(dest .= f.(unpack_args(args, field)...))
@generated copyto_field_res!(dest::Ref, f, args::Tuple, field::Symbol) = :(dest[] = f.(unpack_args(args, field)...))


# Broadcasting implementation
function Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector{Names}}}) where {Names}
    f = bc.f
    args = bc.args
    # Apply broadcast to each field
    res_args = map(Names) do name
        get_field_res(f,args, name)
    end
    HeterogenousVector(NamedTuple{Names}(res_args))
end

@inline Base.@constprop :aggressive function Base.copyto!(dest::HeterogenousVector{T, S},bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector{Names}}, Axes, F, Args}) where {T, S, Names, Axes,F,Args<:Tuple}
    if fieldnames(S) != Names
        throw(ArgumentError("Cannot copy to HeterogenousVector with different field names: $(fieldnames(S)) vs $(Names)"))
    end
    f = bc.f
    args = bc.args
    # Using value types to specialize map_fun is indeed an ugly solution
    # Constant propagation should **usually** make this unnecessary, but
    # benchmarking has shown there are cases where this does not happen (even with aggressive const propagation),
    # causing type unstability and costly runtime dispatch
    function map_fun(::Val{name}) where {name}
        target_field = getfield(NamedTuple(dest), name)
        copyto_field_res!(target_field, f, args, name)
    end
    map(map_fun, Val.(Names))
    dest
end

# Compute segment ranges for each field in the NamedTuple
# The results are zero-indexed ranges, i.e. the first field starts at 0
function _compute_segment_ranges(x::NamedTuple)
    current_idx = 0
    ranges = map(x) do field
        field_length = _field_length(field)
        range = current_idx:(current_idx + field_length - 1)
        current_idx += field_length
        range
    end
    return ranges
end


# Written specifically to deal with cases such as calculate_residuals!() where the destination is an ordinary Array
@inline Base.@constprop :aggressive function Base.copyto!(dest::AbstractArray, bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector{Names}}}) where {Names}
    f = bc.f
    args = bc.args
    hv = find_heterogenous_vector(bc)
    # Points to the first index of the destination array
    dest_idx = firstindex(dest)
    segment_ranges = _compute_segment_ranges(NamedTuple(hv))
    function map_fun(::Val{name}) where {name}
        unpacked_args = unpack_args(args, name)
        segment_range = segment_ranges[name]
        dest_segment = view(dest, dest_idx .+ segment_range)
        dest_segment .= f.(unpacked_args...)
    end
    map(map_fun, Val.(Names))
    return dest
end


# Show methods for HeterogenousVector
Base.summary(hv::HeterogenousVector) = string(typeof(hv), " with members:")
Base.show(io::IO, m::MIME"text/plain", hv::HeterogenousVector) = show(io, m, NamedTuple(hv))

# Copy-catted from RecursiveArrayTools.jl/src/utils.jl
# From Iterators.jl. Moved here since Iterators.jl is not precompile safe anymore.

# Concatenate the output of n iterators
struct Chain{T <: Tuple}
    xss::T
end

# iteratorsize method defined at bottom because of how @generated functions work in 0.6 now

"""
    chain(xs...)

Iterate through any number of iterators in sequence.

```jldoctest
julia> for i in chain(1:3, ['a', 'b', 'c'])
           @show i
       end
i = 1
i = 2
i = 3
i = 'a'
i = 'b'
i = 'c'
```
"""
chain(xss...) = Chain(xss)

Base.length(it::Chain{Tuple{}}) = 0
Base.length(it::Chain) = sum(length, it.xss)

Base.eltype(::Type{Chain{T}}) where {T} = typejoin([eltype(t) for t in T.parameters]...)

function Base.iterate(it::Chain)
    i = 1
    xs_state = nothing
    while i <= length(it.xss)
        xs_state = iterate(it.xss[i])
        xs_state !== nothing && return xs_state[1], (i, xs_state[2])
        i += 1
    end
    return nothing
end

function Base.iterate(it::Chain, state)
    i, xs_state = state
    xs_state = iterate(it.xss[i], xs_state)
    while xs_state == nothing
        i += 1
        i > length(it.xss) && return nothing
        xs_state = iterate(it.xss[i])
    end
    return xs_state[1], (i, xs_state[2])
end

Base.iterate(x::HeterogenousVector) = iterate(Chain(values(NamedTuple(x))))
Base.iterate(x::HeterogenousVector, state) = iterate(Chain(values(NamedTuple(x))), state)
