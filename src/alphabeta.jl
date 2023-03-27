# src/beta.jl

alpha(cstable::BoolCSTable) = 
	CVector(sum(cstable.cs, dims=2)[:], cstable.cindex)

function betasim(T::Type{<:AbstractFloat}, cstable::BoolCSTable)
	csf32 = Float32.(cstable.cs)
	a = alpha(cstable).vector
	b = 1 .- T.(csf32 * csf32') ./ min.(a, a')
	return CMatrix(b, cstable.cindex)
end
betasim(cstable::BoolCSTable) = betasim(Float64, cstable)
