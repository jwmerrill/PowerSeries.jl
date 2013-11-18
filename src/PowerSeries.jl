module PowerSeries

import Base: diff, +, -, *, /, ^, sin, cos, exp, log, show, showcompact

immutable Series{T<:Number} <: Number
  re::T
  ep::Union(T, Series{T})
end

Series(a, b, c...) = Series(a, Series(b, c...))

+(a::Series, b::Series) = Series(a.re + b.re, a.ep + b.ep)
+(a::Number, b::Series) = Series(a + b.re, b.ep)
+(a::Series, b::Number) = Series(a.re + b, a.ep)

-(a::Series) = Series(-a.re, -a.ep)
-(a::Series, b::Series) = Series(a.re - b.re, a.ep - b.ep)
-(a::Number, b::Series) = Series(a - b.re, -b.ep)
-(a::Series, b::Number) = Series(a.re - b, -a.ep)

*{T<:Number}(a::Series{T}, b::Series{T}) = Series(
  a.re*b.re,
  a.re*b.ep + b.re*a.ep +
  restrict(Series(zero(T), a.ep*b.ep))
)
*(a::Number, b::Series) = Series(a*b.re, a*b.ep)
*(a::Series, b::Number) = Series(b*a.re, b*a.ep)

promote_depth{T}(a::Series{T}, b::Series{T}) = Series(b.re, promote_depth(a.ep, b.ep))
promote_depth{T}(a::Series{T}, b::T) = Series(b, promote_depth(a.ep, zero(T)))
promote_depth{T}(a::T, b::Series{T}) = b
promote_depth{T}(a::T, b::T) = b

function /{T<:Number}(a::Series{T}, b::Series{T})
  rb = promote_depth(a.ep, b) # same as restrict(promote_depth(a, b))
  a/b.re - a*pint(diff(b)/(rb*rb))
end
/{T<:Number}(a::Series{T}, b::T) = Series(a.re/b, a.ep/b)
function /{T<:Number}(a::T, b::Series{T})
  rb = restrict(b)
  a/b.re - a*pint(diff(b)/(rb*rb))
end

^{T<:Number}(a::Series{T}, b::Number) = Series(a.re^b, b*a.ep^(b - one(T)))

restrict(p::Series) = _restrict(p.re, p.ep)

_restrict{T<:Number}(r::T, e::T) = r
_restrict{T<:Number}(r::T, e::Series{T}) = Series(r, _restrict(e.re, e.ep))

diff{T<:Number}(p::T) = zero(T)
diff{T<:Number}(p::Series{T}) = _diff(p.ep, one(T))

_diff{T<:Number}(p::T, n::T) = n*p
_diff{T<:Number}(p::Series{T}, n::T) = Series(n*p.re, _diff(p.ep, n + one(T)))

pint{T<:Number}(p::T) = Series(zero(T), p)
pint{T<:Number}(p::Series{T}) = Series(zero(T), _pint(p.re, p.ep, one(T)))

_pint{T<:Number}(r::T, e::T, n::T) = Series(r/n, e/(n + one(T)))
_pint{T<:Number}(r::T, e::Series{T}, n::T) =
  Series(r/n, _pint(e.re, e.ep, n + one(T)))

exp(p::Series) = exp(p.re) + pint(diff(p)*exp(restrict(p)))
sin(p::Series) = sin(p.re) + pint(diff(p)*cos(restrict(p)))
cos(p::Series) = cos(p.re) - pint(diff(p)*sin(restrict(p)))
log(p::Series) = log(p.re) + pint(diff(p)/restrict(p))

function _show_rest{T<:Number}(io::IO, method, r::T, e::Series{T})
  method(r)
  print(io, ",")
  _show_rest(io, method, e.re, e.ep)
end

function _show_rest{T<:Number}(io::IO, method, r::T, e::T)
  method(r)
  print(io, ",")
  method(e)
  print(io, ")")
end

function show{T}(io::IO, p::Series{T})
  print(io, "Series{", T, "}(")
  _show_rest(io, show, p.re, p.ep)
end

function showcompact(io::IO, p::Series)
  print(io, "Series(")
  _show_rest(io, showcompact, p.re, p.ep)
end

export Series, pint, restrict

end