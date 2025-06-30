using Unitful
using Revise
import RecursiveArrayTools


# Copy-catted from DiffEqBase DiffEqBaseUnitfulExt.jl
Value(x::Number) = x
Value(x::Type{T}) where {T <:Number} = T
Value(x::Type{Unitful.AbstractQuantity{T, D, U}}) where {T, D, U} = T
Value(x::Unitful.AbstractQuantity) = x.val


# If we do not specify this as a subtype of AbstractVector, the Broadcast machinery will try to convert it into
# a broadcastable representation, but we make this type to be broadcastable as-is due to our customizations of
# copy and copyto!

# Helper function to wrap scalars in Ref for mutability
_make_mutable(x) = x isa AbstractArray ? x : Ref(x)
_unwrap(x) = x isa Ref ? x[] : x
_set_value!(x::Ref, val) = (x[] = val)
_set_value!(x::AbstractArray, val, idx) = (x[idx] = val)

struct HeterogenousVector{T, S <: NamedTuple} <: AbstractVector{T}
    x::S
    function HeterogenousVector(x::NamedTuple)
        # Wrap scalar fields in Ref for mutability
        mutable_x = map(_make_mutable, x)
        arg_types = map(field -> RecursiveArrayTools.recursive_bottom_eltype(_unwrap(field)), values(mutable_x))
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
        unwrapped_field = _unwrap(field)
        if unwrapped_field isa AbstractArray
            field_length = length(unwrapped_field)
            if current_idx <= idx < current_idx + field_length
                unwrapped_field[idx - current_idx + 1] = val
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

# Update length calculation
Base.length(hv::HeterogenousVector) = sum(field -> begin
    unwrapped = _unwrap(field)
    unwrapped isa AbstractArray ? length(unwrapped) : 1
end, hv.x)

# Complete the HeterogenousVector implementation
Base.length(hv::HeterogenousVector) = sum(field -> field isa AbstractArray ? length(field) : 1, hv.x)
Base.size(hv::HeterogenousVector) = (length(hv),)
Base.firstindex(hv::HeterogenousVector) = 1
Base.lastindex(hv::HeterogenousVector) = length(hv)

function Base.copy(hv::HeterogenousVector)
    copied_x = map(field -> field isa AbstractArray ? copy(field) : _unwrap(field), hv.x)
    HeterogenousVector(copied_x)
end

function Base.copy!(dst::HeterogenousVector, src::HeterogenousVector)
    # Ensure both have the same structure
    if keys(dst.x) != keys(src.x)
        throw(ArgumentError("HeterogenousVectors must have the same field names"))
    end
    
    for name in keys(dst.x)
        src_field = getfield(src.x, name)
        dst_field = getfield(dst.x, name)
        if src_field isa AbstractArray && dst_field isa AbstractArray
            copy!(dst_field, src_field)
        else
            # Update scalar field
            if !(src_field isa Ref) && !(dst_field isa Ref)
                throw(ArgumentError("Fields must be mutable references or arrays"))
            end
            _set_value!(dst_field, _unwrap(src_field))
        end
    end
    return dst
end

function Base.similar(hv::HeterogenousVector{T}) where {T}
    similar_x = map(field -> field isa AbstractArray ? similar(field) : zero(typeof(_unwrap(field))), hv.x)
    HeterogenousVector(similar_x)
end

function Base.similar(hv::HeterogenousVector, ::Type{S}) where {S}
    similar_x = map(field -> field isa AbstractArray ? similar(field, S) : zero(_unwrap(S)), hv.x)
    HeterogenousVector(similar_x)
end

function Base.zero(hv::HeterogenousVector)
    zero_x = map(field -> field isa AbstractArray ? zero(field) : zero(_unwrap(field)), hv.x)
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

# Broadcasting implementation
function Base.copy(bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector}})
    hv = find_heterogenous_vector(bc)
    f = bc.f
    args = bc.args
    
    # Apply broadcast to each field
    result_x = map(keys(hv.x)) do name
        field_args = map(arg -> arg isa HeterogenousVector ? getfield(arg.x, name) : arg, args)
        f.(field_args...)
    end
    
    HeterogenousVector(NamedTuple{keys(hv.x)}(result_x))
end

function Base.copyto!(dest::HeterogenousVector, bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector}})
    f = bc.f
    args = bc.args
    
    for name in keys(dest.x)
        field_args = map(arg -> arg isa HeterogenousVector ? getfield(arg.x, name) : arg, args)
        dest_field = getfield(dest.x, name)
        
        if dest_field isa AbstractArray
            dest_field .= f.(field_args...)
        else
            # Update scalar field
            new_value = f(field_args...)
            # Ensure the field is mutable
            if !(dest_field isa Ref)
                throw(ArgumentError("Fields must be mutable references or arrays"))
            end
            # Set the new value
            _set_value!(dst_field, new_value)
        end
    end
    
    return dest
end

# Show methods for HeterogenousVector
Base.summary(hv::HeterogenousVector) = string(typeof(hv), " with members:")
Base.show(io::IO, m::MIME"text/plain", hv::HeterogenousVector) = show(io, m, hv.x)
