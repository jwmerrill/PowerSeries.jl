module PowerSeries

import Base: diff, +, -, *, /, sin, cos, exp

immutable Series{T<:Real} <: Number
  re::T
  ep::Union(T, Series{T})
end

+(a::Series, b::Series) = Series(a.re + b.re, a.ep + b.ep)
+(a::Real, b::Series) = Series(a + b.re, b.ep)
+(a::Series, b::Real) = Series(a.re + b, a.ep)

-(a::Series) = Series(-a.re, -a.ep)
-(a::Series, b::Series) = Series(a.re - b.re, a.ep - b.ep)
-(a::Real, b::Series) = Series(a - b.re, -b.ep)
-(a::Series, b::Real) = Series(a.re - b, -a.ep)

*{T<:Real}(a::Series{T}, b::Series{T}) = Series(
  a.re*b.re,
  a.re*b.ep + b.re*a.ep +
  restrict(Series(zero(T), a.ep*b.ep))
)
*(a::Real, b::Series) = Series(a*b.re, a*b.ep)
*(a::Series, b::Real) = Series(b*a.re, b*a.ep)

restrict(p::Series) = _restrict(p.re, p.ep)

_restrict{T<:Real}(r::T, e::T) = r
_restrict{T<:Real}(r::T, e::Series{T}) = Series(r, _restrict(e.re, e.ep))

diff{T<:Real}(p::T) = zero(T)
diff{T<:Real}(p::Series{T}) = _diff(p.ep, one(T))

_diff{T<:Real}(p::T, n::T) = n*p
_diff{T<:Real}(p::Series{T}, n::T) = Series(n*p.re, _diff(p.ep, n + one(T)))

pint{T<:Real}(p::T) = Series(zero(T), p)
pint{T<:Real}(p::Series{T}) = Series(zero(T), _pint(p.re, p.ep, one(T)))

_pint{T<:Real}(r::T, e::T, n::T) = Series(r/n, e/(n + one(T)))
_pint{T<:Real}(r::T, e::Series{T}, n::T) =
  Series(r/n, _pint(e.re, e.ep, n + one(T)))

exp(p::Series) = exp(p.re) + pint(diff(p)*exp(restrict(p)))
sin(p::Series) = sin(p.re) + pint(diff(p)*cos(restrict(p)))
cos(p::Series) = cos(p.re) - pint(diff(p)*sin(restrict(p)))

export Series, pint, restrict

end