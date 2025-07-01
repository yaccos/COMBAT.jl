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


# Custom property access for clean external interface
function Base.getproperty(hv::HeterogenousVector, name::Symbol)
    if name === :x
        # Allow access to the internal NamedTuple
        return getfield(hv, :x)
    elseif haskey(hv.x, name)
        # Access field and unwrap if it's a Ref
        field = getfield(hv.x, name)
        return _unwrap(field)
    else
        error("HeterogenousVector has no field $name")
    end
end

function Base.setproperty!(hv::HeterogenousVector, name::Symbol, value)
    if name === :x
        # Prevent direct assignment to internal structure
        error("Cannot directly assign to internal field :x")
    elseif haskey(hv.x, name)
        # Set field value, wrapping in Ref if it's a scalar
        field = getfield(hv.x, name)
        _set_value!(field, value)
    else
        error("HeterogenousVector has no field $name")
    end
end

function Base.propertynames(hv::HeterogenousVector)
    return keys(hv.x)
end






function Base.getindex(hv::HeterogenousVector{T}, idx::Int) where {T}
    current_idx = 1
    for (name, field) in pairs(hv.x)
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

function Base.setindex!(hv::HeterogenousVector{T}, val, idx::Int) where {T}
    current_idx = 1
    for (name, field) in pairs(hv.x)
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
Base.length(hv::HeterogenousVector) = sum(_field_length, hv.x)

Base.size(hv::HeterogenousVector) = (length(hv),)
Base.firstindex(hv::HeterogenousVector) = 1
Base.lastindex(hv::HeterogenousVector) = length(hv)

_copy_field(field::Ref) = _unwrap(field) |> Ref
_copy_field(field::AbstractArray) = copy(field)
function Base.copy(hv::HeterogenousVector)
    copied_x = map(_copy_field, hv.x)
    HeterogenousVector(copied_x)
end



_copy_field!(dst::Ref, src::Ref) = _set_value!(dst, _unwrap(src))
_copy_field!(dst::AbstractArray, src::AbstractArray) = copy!(dst, src)
# Update the copy! function to handle the new interface properly
function Base.copy!(dst::HeterogenousVector, src::HeterogenousVector)
    # Ensure both have the same structure
    if keys(dst.x) != keys(src.x)
        throw(ArgumentError("HeterogenousVectors must have the same field names"))
    end
    
    for name in keys(dst.x)
        src_field = getfield(src.x, name)
        dst_field = getfield(dst.x, name)
        
        _copy_field!(dst_field, src_field)
    end
    return dst
end

Base.zero(x::Ref) = _unwrap(x) |> zero
Base.similar(x::Ref) = zero(x)

function Base.similar(hv::HeterogenousVector{T}) where {T}
    similar_x = map(similar, hv.x)
    HeterogenousVector(similar_x)
end

function Base.zero(hv::HeterogenousVector)
    zero_x = map(zero, hv.x)
    HeterogenousVector(zero_x)
end

# Broadcasting support for HeterogenousVector
Base.BroadcastStyle(::Type{<:HeterogenousVector}) = Broadcast.Style{HeterogenousVector}()
Base.BroadcastStyle(::Broadcast.Style{HeterogenousVector}, ::Base.Broadcast.BroadcastStyle) = Broadcast.Style{HeterogenousVector}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector}})
    hv = find_heterogenous_vector(bc)
    similar(hv)
end

# Helper function to find HeterogenousVector in broadcast arguments
find_heterogenous_vector(bc::Base.Broadcast.Broadcasted) = find_heterogenous_vector(bc.args)
find_heterogenous_vector(args::Tuple) = find_heterogenous_vector(find_heterogenous_vector(args[1]), Base.tail(args))
find_heterogenous_vector(x) = x
find_heterogenous_vector(::Tuple{}) = nothing
find_heterogenous_vector(x::HeterogenousVector, rest) = x
find_heterogenous_vector(::Any, rest) = find_heterogenous_vector(rest)

# Generic field unpacking for HeterogenousVector
@generated function unpack_field(hv::HeterogenousVector{T, S}, field::Symbol) where {T, S}
    if field in fieldnames(S)
        :(getfield(hv.x, field))
    else
        :(throw(ArgumentError("Field $field not found in HeterogenousVector")))
    end
end



_get_field_arg(arg::HeterogenousVector, name::Symbol) = getfield(arg.x, name) |> _unwrap
_get_field_arg(arg::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector}}, name::Symbol) = Broadcast.broadcasted(arg.f,_get_field_arg(arg.args,name)...)
# Fix the tuple handling for _get_field_arg
_get_field_arg(args::Tuple, name::Symbol) = map(arg -> _get_field_arg(arg, name), args)

_get_field_arg(arg, name::Symbol) = arg
_copyto_dest_field!(dest_field::AbstractArray,field_args, f) = (dest_field .= f.(field_args...))
_copyto_dest_field!(dest_field::Ref,field_args, f) = _set_value!(dest_field,f(field_args...))
# If the field is an array, apply the function element-wise
_get_field_res(field_args, f, ::AbstractArray) = f.(field_args...)
# If the field is a scalar, apply the function directly
_get_field_res(field_args, f, ::Any) = f(field_args...)
# Update broadcasting copyto! to use the variable name correctly

function Base.copyto!(dest::HeterogenousVector, bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector}})
    f = bc.f
    args = bc.args
    
    for name in keys(dest.x)
        field_args = map(arg -> _get_field_arg(arg,name), args)
        dest_field = getfield(dest.x, name)
        _copyto_dest_field!(dest_field, field_args, f)
    end
    return dest
end

# Broadcasting implementation
function Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector}})
    hv = find_heterogenous_vector(bc)
    f = bc.f
    args = bc.args
    
    # Apply broadcast to each field
    result_x = map(keys(hv.x)) do name
        field_args = map(arg -> _get_field_arg(arg, name), args)
        proto_dest_field = getfield(hv.x, name)
        _get_field_res(field_args, f, proto_dest_field)
    end
    HeterogenousVector(NamedTuple{keys(hv.x)}(result_x))
end



# Show methods for HeterogenousVector
Base.summary(hv::HeterogenousVector) = string(typeof(hv), " with members:")
Base.show(io::IO, m::MIME"text/plain", hv::HeterogenousVector) = show(io, m, hv.x)
