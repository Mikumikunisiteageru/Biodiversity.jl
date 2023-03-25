# src/beta.jl

alpha(cstable::DiscreteCSTable) = 
	CVector(sum(cstable.cs, dims=2)[:], cstable.cindex)

function betasim(T::Type{<:AbstractFloat}, cstable::DiscreteCSTable)
	csf32 = Float32.(cstable.cs)
	a = alpha(cstable).vector
	b = 1 .- T.(csf32 * csf32') ./ min.(a, a')
	return CMatrix(b, cstable.cindex)
end
betasim(cstable::DiscreteCSTable) = betasim(Float64, cstable)
