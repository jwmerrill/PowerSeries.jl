module PowerSeries

import Base: diff, +, -, *, /, sin, cos, exp

immutable PowerSeriesType{T<:Real} <: Number
  re::T
  ep::Union(T, PowerSeriesType{T})
end

+(a::PowerSeriesType, b::PowerSeriesType) = PowerSeriesType(a.re + b.re, a.ep + b.ep)
+(a::Real, b::PowerSeriesType) = PowerSeriesType(a + b.re, b.ep)
+(a::PowerSeriesType, b::Real) = PowerSeriesType(a.re + b, a.ep)

-(a::PowerSeriesType) = PowerSeriesType(-a.re, -a.ep)
-(a::PowerSeriesType, b::PowerSeriesType) = PowerSeriesType(a.re - b.re, a.ep - b.ep)
-(a::Real, b::PowerSeriesType) = PowerSeriesType(a - b.re, b.ep)
-(a::PowerSeriesType, b::Real) = PowerSeriesType(a.re - b, a.ep)

*{T<:Real}(a::PowerSeriesType{T}, b::PowerSeriesType{T}) = PowerSeriesType(
  a.re*b.re,
  a.re*b.ep + b.re*a.ep +
  restrict(PowerSeriesType(zero(T), a.ep*b.ep))
)
*(a::Real, b::PowerSeriesType) = PowerSeriesType(a*b.re, a*b.ep)
*(a::PowerSeriesType, b::Real) = PowerSeriesType(b*a.re, b*a.ep)

restrict(p::PowerSeriesType) = _restrict(p.re, p.ep)

_restrict{T<:Real}(r::T, e::T) = r
_restrict{T<:Real}(r::T, e::PowerSeriesType{T}) = PowerSeriesType(r, _restrict(e.re, e.ep))

diff{T<:Real}(p::T) = zero(T)

diff{T<:Real}(p::PowerSeriesType{T}) = _diff(p.ep, one(T))

function _diff{T<:Real}(p::T, n::T)
  n*p
end

function _diff{T<:Real}(p::PowerSeriesType{T}, n::T)
  PowerSeriesType(n*p.re, _diff(p.ep, n + one(T)))
end

pint{T<:Real}(p::T) = PowerSeriesType(zero(T), p)

pint{T<:Real}(p::PowerSeriesType{T}) = PowerSeriesType(zero(T), _pint(p.re, p.ep, one(T)))

function _pint{T<:Real}(r::T, e::T, n::T)
  PowerSeriesType(r/n, e/(n + one(T)))
end

function _pint{T<:Real}(r::T, e::PowerSeriesType{T}, n::T)
  PowerSeriesType(r/n, _pint(e.re, e.ep, n + one(T)))
end

exp(p::PowerSeriesType) = exp(p.re) + pint(diff(p)*exp(restrict(p)))
sin(p::PowerSeriesType) = sin(p.re) + pint(diff(p)*cos(restrict(p)))
cos(p::PowerSeriesType) = cos(p.re) - pint(diff(p)*sin(restrict(p)))

export PowerSeriesType, pint, restrict

end