module PowerSeries

import Base: sin, cos, exp, log

abstract AbstractSeries{T<:Real, N} <: Number

immutable Series1{T} <: AbstractSeries{T, 1}
  c0::T
  c1::T
end

immutable Series2{T} <: AbstractSeries{T, 2}
  c0::T
  c1::T
  c2::T
end

immutable Series3{T} <: AbstractSeries{T, 3}
  c0::T
  c1::T
  c2::T
  c3::T
end

series_types = [Series1, Series2, Series3]

function series(n...)
  l = length(n)
  series_types[l - 1](promote(n...)...)
end

restrict(x::Series1) = x.c0
restrict(x::Series2) = Series1(x.c0, x.c1)
restrict(x::Series3) = Series2(x.c0, x.c1, x.c2)

constant(x::Series1) = x.c0
constant(x::Series2) = x.c0
constant(x::Series3) = x.c0
constant(x::Real) = x

# Horner's method
polyval(x::Series1, eps::Real) = x.c0 + eps*x.c1
polyval(x::Series2, eps::Real) = x.c0 + eps*(x.c1 + eps*x.c2)
polyval(x::Series3, eps::Real) = x.c0 + eps*(x.c1 + eps*(x.c2 + eps*x.c3))

polyder(x::Real) = zero(typeof(x))
polyder(x::Series1) = x.c1
polyder(x::Series2) = Series1(x.c1, 2*x.c2)
polyder(x::Series3) = Series2(x.c1, 2*x.c2, 3*x.c3)

polyint(x::Real) = Series1(zero(typeof(x)), x)
polyint{T}(x::Series1{T}) = Series2(zero(T), x.c0, x.c1/2)
polyint{T}(x::Series2{T}) = Series3(zero(T), x.c0, x.c1/2, x.c2/3)

+(x::Series1, y::Series1) = Series1(x.c0 + y.c0, x.c1 + y.c1)
+(x::Series2, y::Series2) = Series2(x.c0 + y.c0, x.c1 + y.c1, x.c2 + y.y2)
+(x::Series3, y::Series3) = Series3(x.c0 + y.c0, x.c1 + y.c1, x.c2 + y.c2, x.c3 + y.c3)

+{T<:Real, S}(c::T, x::Series1{S}) = Series1{promote_type(T, S)}(c + x.c0, x.c1)
+{T<:Real, S}(c::T, x::Series2{S}) = Series2{promote_type(T, S)}(c + x.c0, x.c1, x.c2)
+{T<:Real, S}(c::T, x::Series3{S}) = Series3{promote_type(T, S)}(c + x.c0, x.c1, x.c2, x.c3)

+{T<:Real, S}(x::Series1{S}, c::T) = Series1{promote_type(T, S)}(x.c0 + c, x.c1)
+{T<:Real, S}(x::Series2{S}, c::T) = Series2{promote_type(T, S)}(x.c0 + c, x.c1, x.c2)
+{T<:Real, S}(x::Series3{S}, c::T) = Series3{promote_type(T, S)}(x.c0 + c, x.c1, x.c2, x.c3)

-(x::Series1) = Series1(-x.c0, -x.c1)
-(x::Series2) = Series2(-x.c0, -x.c1, -x.c2)
-(x::Series3) = Series3(-x.c0, -x.c1, -x.c2, -x.c3)

-(x::Series1, y::Series1) = Series1(x.c0 - y.c0, x.c1 - y.c1)
-(x::Series2, y::Series2) = Series2(x.c0 - y.c0, x.c1 - y.c1, x.c2 - y.y2)
-(x::Series3, y::Series3) = Series3(x.c0 - y.c0, x.c1 - y.c1, x.c2 - y.c2, x.c3 - y.c3)

-{T<:Real, S}(c::T, x::Series1{S}) = Series1{promote_type(T, S)}(c - x.c0, -x.c1)
-{T<:Real, S}(c::T, x::Series2{S}) = Series2{promote_type(T, S)}(c - x.c0, -x.c1, -x.c2)
-{T<:Real, S}(c::T, x::Series3{S}) = Series3{promote_type(T, S)}(c - x.c0, -x.c1, -x.c2, -x.c3)

-{T<:Real, S}(x::Series1{S}, c::T) = Series1{promote_type(T, S)}(x.c0 - c, x.c1)
-{T<:Real, S}(x::Series2{S}, c::T) = Series2{promote_type(T, S)}(x.c0 - c, x.c1, x.c2)
-{T<:Real, S}(x::Series3{S}, c::T) = Series3{promote_type(T, S)}(x.c0 - c, x.c1, x.c2, x.c3)

*(x::Series1, y::Series1) = Series1(x.c0*y.c0, x.c0*y.c1 + x.c1*y.c0)
*(x::Series2, y::Series2) = Series2(x.c0*y.c0, x.c0*y.c1 + x.c1*y.c0, x.c0*y.c2 + x.c1*y.c1 + x.c2*y.c0)
*(x::Series3, y::Series3) = Series3(
  x.c0*y.c0,
  x.c0*y.c1 + x.c1*y.c0,
  x.c0*y.c2 + x.c1*y.c1 + x.c2*y.c0,
  x.c0*y.c3 + x.c1*y.c2 + x.c2*y.c1 + x.c3*y.c0
)

*(c::Real, x::Series1) = Series1(c*x.c0, c*x.c1)
*(c::Real, x::Series2) = Series2(c*x.c0, c*x.c1, c*x.c2)
*(c::Real, x::Series3) = Series3(c*x.c0, c*x.c1, c*x.c2, c*x.c3)

*(x::Series1, c::Real) = Series1(x.c0*c, x.c1*c)
*(x::Series2, c::Real) = Series2(x.c0*c, x.c1*c, x.c2*c)
*(x::Series3, c::Real) = Series3(x.c0*c, x.c1*c, x.c2*c, x.c3*c)

function /{T, S, N}(x::AbstractSeries{T, N}, y::AbstractSeries{S, N})
  ry = restrict(y)
  constant(x)/constant(y) + polyint(polyder(x)/ry - restrict(x)*polyder(y)/(ry*ry))
end

function /(c::Real, y::AbstractSeries)
  ry = restrict(y)
  c/constant(y) - polyint(c*polyder(y)/(ry*ry))
end

/(x::Series1, c::Real) = Series1(x.c0/c, x.c1/c)
/(x::Series2, c::Real) = Series2(x.c0/c, x.c1/c, x.c2/c)
/(x::Series3, c::Real) = Series3(x.c0/c, x.c1/c, x.c2/c, x.c3/c)

function ^{T, S, N}(x::AbstractSeries{T, N}, y::AbstractSeries{S, N})
  rx = restrict(x)
  ry = restrict(y)
  constant(x)^constant(y) + polyint(polyder(x)*rx^(ry - 1)*ry + log(rx)*rx^ry*polyder(y))
end

_series_pow_const(x, y) = constant(x)^y + polyint(polyder(x)*restrict(x)^(y - 1)*y)

# First two are to fix redundancy warnings
^(x::AbstractSeries, y::Rational) = _series_pow_const(x, y)
^(x::AbstractSeries, y::Integer) = _series_pow_const(x, y)
^(x::AbstractSeries, y::Real) = _series_pow_const(x, y)

^(::MathConst{:e}, y::AbstractSeries) = exp(y)
^(x::Real, y::AbstractSeries) = x^constant(y) + polyint(log(x)*x^restrict(y)*polyder(y))

sin(x::AbstractSeries) = sin(constant(x)) + polyint(polyder(x)*cos(restrict(x)))
cos(x::AbstractSeries) = cos(constant(x)) - polyint(polyder(x)*sin(restrict(x)))
exp(x::AbstractSeries) = exp(constant(x)) + polyint(polyder(x)*exp(restrict(x)))
log(x::AbstractSeries) = log(constant(x)) + polyint(polyder(x)/restrict(x))

export series, restrict, constant, polyint, polyval, polyder

end