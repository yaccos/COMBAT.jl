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
mutable struct HeterogenousVector{T, S <: NamedTuple} <: AbstractVector{T}
    data::S
    function HeterogenousVector(data::NamedTuple)
        arg_types = map(RecursiveArrayTools.recursive_bottom_eltype,values(data))
        T = promote_type(arg_types...)
        new{T,typeof(data)}(data)
    end

    # Constructor with keyword arguments - field names determined by keywords
    function HeterogenousVector(; kwargs...)
        data = NamedTuple(kwargs)
        arg_types = map(RecursiveArrayTools.recursive_bottom_eltype, values(data))
        T = promote_type(arg_types...)
        new{T, typeof(data)}(data)
    end

    # Constructor with positional arguments (keeping the original functionality)
    function HeterogenousVector(args...)
        # Convert positional args to NamedTuple with generated names
        names = ntuple(i -> Symbol("field_$i"), length(args))
        data = NamedTuple{names}(args)
        HeterogenousVector(data)
    end
end

# Complete the HeterogenousVector implementation
Base.length(hv::HeterogenousVector) = sum(field -> field isa AbstractArray ? length(field) : 1, hv.data)
Base.size(hv::HeterogenousVector) = (length(hv),)
Base.firstindex(hv::HeterogenousVector) = 1
Base.lastindex(hv::HeterogenousVector) = length(hv)

function Base.getindex(hv::HeterogenousVector{T}, idx::Int) where {T}
    current_idx = 1
    for (name, field) in pairs(hv.data)
        if field isa AbstractArray
            field_length = length(field)
            if current_idx <= idx < current_idx + field_length
                return field[idx - current_idx + 1]
            end
            current_idx += field_length
        else
            if idx == current_idx
                return field
            end
            current_idx += 1
        end
    end
    throw(BoundsError(hv, idx))
end

function Base.setindex!(hv::HeterogenousVector{T}, val, idx::Int) where {T}
    current_idx = 1
    for (name, field) in pairs(hv.data)
        if field isa AbstractArray
            field_length = length(field)
            if current_idx <= idx < current_idx + field_length
                field[idx - current_idx + 1] = val
                return val
            end
            current_idx += field_length
        else
            if idx == current_idx
                # For scalar fields, we need to update the NamedTuple
                new_data = merge(hv.data, NamedTuple{(name,)}((val,)))
                hv.data = new_data
                return val
            end
            current_idx += 1
        end
    end
    throw(BoundsError(hv, idx))
end

function Base.copy(hv::HeterogenousVector)
    copied_data = map(field -> field isa AbstractArray ? copy(field) : field, hv.data)
    HeterogenousVector(copied_data)
end

function Base.copy!(dst::HeterogenousVector, src::HeterogenousVector)
    # Ensure both have the same structure
    if keys(dst.data) != keys(src.data)
        throw(ArgumentError("HeterogenousVectors must have the same field names"))
    end
    
    for name in keys(dst.data)
        src_field = getfield(src.data, name)
        dst_field = getfield(dst.data, name)
        if src_field isa AbstractArray && dst_field isa AbstractArray
            copy!(dst_field, src_field)
        else
            # Update scalar field
            new_data = merge(dst.data, NamedTuple{(name,)}((src_field,)))
            dst.data = new_data
        end
    end
    return dst
end

function Base.similar(hv::HeterogenousVector{T}) where {T}
    similar_data = map(field -> field isa AbstractArray ? similar(field) : zero(typeof(field)), hv.data)
    HeterogenousVector(similar_data)
end

function Base.similar(hv::HeterogenousVector, ::Type{S}) where {S}
    similar_data = map(field -> field isa AbstractArray ? similar(field, S) : zero(S), hv.data)
    HeterogenousVector(similar_data)
end

function Base.zero(hv::HeterogenousVector)
    zero_data = map(field -> field isa AbstractArray ? zero(field) : zero(field), hv.data)
    HeterogenousVector(zero_data)
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
        :(getfield(hv.data, field))
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
    result_data = map(keys(hv.data)) do name
        field_args = map(arg -> arg isa HeterogenousVector ? getfield(arg.data, name) : arg, args)
        f.(field_args...)
    end
    
    HeterogenousVector(NamedTuple{keys(hv.data)}(result_data))
end

function Base.copyto!(dest::HeterogenousVector, bc::Broadcast.Broadcasted{Broadcast.Style{HeterogenousVector}})
    f = bc.f
    args = bc.args
    
    for name in keys(dest.data)
        field_args = map(arg -> arg isa HeterogenousVector ? getfield(arg.data, name) : arg, args)
        dest_field = getfield(dest.data, name)
        
        if dest_field isa AbstractArray
            dest_field .= f.(field_args...)
        else
            # Update scalar field
            new_value = f(field_args...)
            new_data = merge(dest.data, NamedTuple{(name,)}((new_value,)))
            dest.data = new_data
        end
    end
    
    return dest
end

# Show methods for HeterogenousVector
function Base.show(io::IO, hv::HeterogenousVector)
    print(io, "HeterogenousVector with fields: ")
    print(io, join(string.(keys(hv.data)), ", "))
end

function Base.show(io::IO, mime::MIME"text/plain", hv::HeterogenousVector)
    multiline = get(io, :multiline, true)
    if multiline
        println(io, "HeterogenousVector with $(length(hv.data)) fields and $(length(hv)) elements:")
        for (name, field) in pairs(hv.data)
            print(io, "  $name: ")
            if field isa AbstractArray
                println(io, "$(typeof(field)) with $(length(field)) elements:")
                for (i, value) in enumerate(field)
                    println(io, "    [$i] = $value")
                end
            else
                println(io, "$field")
            end
        end
    else
        Base.show_default(io, hv)
    end
end