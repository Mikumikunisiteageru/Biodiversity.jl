# src/types.jl

struct StrIndex
	list::Vector{String}
	dict::Dict{String, Int}
	len::Int
	function StrIndex(list::AbstractVector{<:AbstractString})
		allunique(list) || 
			throw(ArgumentError("some string elements have duplicates"))
		len = length(list)
		dict = Dict(list .=> 1:len)
		return new(list, dict, len)
	end
end

Base.length(strindex::StrIndex) = strindex.len

Base.getindex(strindex::StrIndex, i::Integer) = strindex.list[i]
Base.getindex(strindex::StrIndex, str::AbstractString) = strindex.dict[str]

Base.haskey(strindex::StrIndex, key::AbstractString) = 
	haskey(strindex.dict, key)

struct CVector{T<:Real} <: AbstractVector{T}
	vector::Vector{T}
	cindex::StrIndex
	function CVector(vector::AbstractVector{T}, cindex::StrIndex) where 
			{T <: Real}
		length(vector) == length(cindex) || 
			throw(DimensionMismatch("the numbers of communities do not match"))
		return new{T}(vector, cindex)
	end
end
Base.length(v::CVector) = length(v.vector)
Base.size(v::CVector) = size(v.vector)
Base.getindex(v::CVector, i...) = getindex(v.vector, i...)

struct CMatrix{T<:Real} <: AbstractMatrix{T}
	matrix::Matrix{T}
	cindex::StrIndex
	function CMatrix(matrix::AbstractMatrix{T}, cindex::StrIndex) where 
			{T <: Real}
		size(matrix, 1) == size(matrix, 2) == length(cindex) || 
			throw(DimensionMismatch("the numbers of communities do not match"))
		return new{T}(matrix, cindex)
	end
end
Base.size(m::CMatrix) = size(m.matrix)
Base.getindex(m::CMatrix, i...) = getindex(m.matrix, i...)

abstract type AbstractCSTable end

function check_s_indexed(cstable::AbstractCSTable, s::AbstractString)
	haskey(cstable.sindex, s) || 
		throw(KeyError("the species $s not indexed in the table"))
	return
end

function check_c_indexed(cstable::AbstractCSTable, c::AbstractString)
	haskey(cstable.cindex, c) || 
		throw(KeyError("the community $c not indexed in the table"))
	return
end

struct DiscreteCSTable <: AbstractCSTable
	cs::Matrix{Float32}
	cindex::StrIndex
	sindex::StrIndex
	function DiscreteCSTable(cs::AbstractMatrix{<:Real}, 
			cindex::StrIndex, sindex::StrIndex)
		size(cs, 1) == length(cindex) || 
			throw(DimensionMismatch("the numbers of communities do not match"))
		size(cs, 2) == length(sindex) || 
			throw(DimensionMismatch("the numbers of species do not match"))
		csf32 = Float32.(Bool.(cs)) # may throw an InexactError
		return new(csf32, cindex, sindex)
	end
end

function DiscreteCSTable(cs::AbstractMatrix{<:Real}, 
			cindex::AbstractVector{<:AbstractString}, 
			sindex::AbstractVector{<:AbstractString})
	return DiscreteCSTable(cs, StrIndex(cindex), StrIndex(sindex))
end
	
function DiscreteCSTable(csrecords::AbstractVector{NTuple{2, String}})
	cindex = StrIndex(sort!(unique!(first.(csrecords))))
	sindex = StrIndex(sort!(unique!(last.(csrecords))))
	cs = zeros(Float32, cindex.len, sindex.len)
	for (c, s) = csrecords
		cs[cindex[c], sindex[s]] = 1f0
	end
	return DiscreteCSTable(cs, cindex, sindex)
end

function DiscreteCSTable(filename::AbstractString; 
		delim::AbstractChar='\t', sfirst::Bool=false)
	rectable = readdlm(filename, delim, String)
	size(rectable, 2) == 2 || throw(
		DimensionMismatch("the TSV file does not contain exactly two columns"))
	return sfirst ? 
		DiscreteCSTable(reverse.(Tuple.(eachrow(rectable)))) : 
		DiscreteCSTable(Tuple.(eachrow(rectable)))
end	

function ispresent(T::Type{<:Real}, cstable::DiscreteCSTable, 
		s::AbstractString, c::AbstractString)
	check_s_indexed(cstable, s)
	check_c_indexed(cstable, c)
	return T(cstable.cs[cstable.cindex[c], cstable.sindex[s]])
end
ispresent(cstable::DiscreteCSTable, s::AbstractString, c::AbstractString) = 
	ispresent(Bool, cstable, s, c)

function findssfromc(cstable::DiscreteCSTable, c::AbstractString; 
		present::Bool=true)
	check_c_indexed(cstable, c)
	return cstable.sindex.list[cstable.cs[cstable.cindex[c], :] .== present]
end

function findccfroms(cstable::DiscreteCSTable, s::AbstractString; 
		present::Bool=true)
	check_s_indexed(cstable, s)
	return cstable.cindex.list[:, cstable.cs[cstable.sindex[s]] .== present]
end
