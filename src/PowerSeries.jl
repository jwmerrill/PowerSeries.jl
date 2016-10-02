module PowerSeries

import Base: convert, sqrt, exp, log, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, csc, sec, cot, acsc, asec, acot, csch, sech, coth, acsch, asech, acoth, gamma, polygamma, floor, ceil, round, sign, abs, <, /, *, ^, +, -

abstract AbstractSeries{T<:Real, N} <: Number

include("generate_type.jl")

series_types = DataType[]

function generate(order)
  while order > length(series_types)
    push!(series_types, generate_type(length(series_types) + 1))
  end
end

constant(x::Real) = x
polyder(x::Real) = zero(typeof(x))

function /{T, S, N}(x::AbstractSeries{T, N}, y::AbstractSeries{S, N})
  ry = restrict(y)
  constant(x)/constant(y) + polyint(polyder(x)/ry - restrict(x)*polyder(y)/(ry*ry))
end

function /(c::Real, y::AbstractSeries)
  ry = restrict(y)
  c/constant(y) - polyint(c*polyder(y)/(ry*ry))
end

function ^{T, S, N}(x::AbstractSeries{T, N}, y::AbstractSeries{S, N})
  rx = restrict(x)
  ry = restrict(y)
  constant(x)^constant(y) + polyint(polyder(x)*rx^(ry - 1)*ry + log(rx)*rx^ry*polyder(y))
end

function _series_pow_const(x, y)
  if y == 0
    return one(x)
  else
    return constant(x)^y + polyint(polyder(x)*restrict(x)^(y - 1)*y)
  end
end

# First two are to fix redundancy warnings
^(x::AbstractSeries, y::Rational) = _series_pow_const(x, y)
^(x::AbstractSeries, y::Integer) = _series_pow_const(x, y)
^(x::AbstractSeries, y::Real) = _series_pow_const(x, y)

^(::Irrational{:e}, y::AbstractSeries) = exp(y)
^(x::Real, y::AbstractSeries) = x^constant(y) + polyint(log(x)*x^restrict(y)*polyder(y))

sqrt(x::AbstractSeries) = sqrt(constant(x)) + polyint(polyder(x)/(2*sqrt(restrict(x))))
exp(x::AbstractSeries) = exp(constant(x)) + polyint(polyder(x)*exp(restrict(x)))
log(x::AbstractSeries) = log(constant(x)) + polyint(polyder(x)/restrict(x))
sin(x::AbstractSeries) = sin(constant(x)) + polyint(polyder(x)*cos(restrict(x)))
cos(x::AbstractSeries) = cos(constant(x)) - polyint(polyder(x)*sin(restrict(x)))
tan(x::AbstractSeries) = tan(constant(x)) + polyint(polyder(x)*sec(restrict(x))^2)
asin(x::AbstractSeries) = asin(constant(x)) + polyint(polyder(x)/sqrt(1 - restrict(x)^2))
acos(x::AbstractSeries) = acos(constant(x)) - polyint(polyder(x)/sqrt(1 - restrict(x)^2))
atan(x::AbstractSeries) = atan(constant(x)) + polyint(polyder(x)/(1 + restrict(x)^2))
sinh(x::AbstractSeries) = sinh(constant(x)) + polyint(polyder(x)*cosh(restrict(x)))
cosh(x::AbstractSeries) = cosh(constant(x)) + polyint(polyder(x)*sinh(restrict(x)))
tanh(x::AbstractSeries) = tanh(constant(x)) + polyint(polyder(x)*sech(restrict(x))^2)
asinh(x::AbstractSeries) = asinh(constant(x)) + polyint(polyder(x)/sqrt(restrict(x)^2 + 1))
acosh(x::AbstractSeries) = acosh(constant(x)) + polyint(polyder(x)/sqrt(restrict(x)^2 - 1))
atanh(x::AbstractSeries) = atanh(constant(x)) + polyint(polyder(x)/(1 - restrict(x)^2))
csc(x::AbstractSeries) = csc(constant(x)) - polyint(polyder(x)*cot(restrict(x))*csc(restrict(x)))
sec(x::AbstractSeries) = sec(constant(x)) + polyint(polyder(x)*tan(restrict(x))*sec(restrict(x)))
cot(x::AbstractSeries) = cot(constant(x)) - polyint(polyder(x)*csc(restrict(x))^2)

acsc(x::AbstractSeries) =
  acsc(constant(x)) - polyint(polyder(x)/(restrict(x)^2*sqrt(1-1/restrict(x)^2)))

asec(x::AbstractSeries) =
  asec(constant(x)) + polyint(polyder(x)/(restrict(x)^2*sqrt(1-1/restrict(x)^2)))

acot(x::AbstractSeries) = acot(constant(x)) - polyint(polyder(x)/(1 + restrict(x)^2))

csch(x::AbstractSeries) =
  csch(constant(x)) - polyint(polyder(x)*coth(restrict(x))*csch(restrict(x)))

sech(x::AbstractSeries) =
  sech(constant(x)) - polyint(polyder(x)*tanh(restrict(x))*sech(restrict(x)))

coth(x::AbstractSeries) = coth(constant(x)) - polyint(polyder(x)*csch(restrict(x))^2)

acsch(x::AbstractSeries) =
  acsch(constant(x)) - polyint(polyder(x)/(restrict(x)^2*sqrt(1+1/restrict(x)^2)))

function asech(x::AbstractSeries)
  rx = restrict(x)
  asech(constant(x)) + polyint(polyder(x)*sqrt((1 - rx)/(1 + rx))/(rx*(rx - 1)))
end

acoth(x::AbstractSeries) = acoth(constant(x)) + polyint(polyder(x)/(1 - restrict(x)^2))

gamma(x::AbstractSeries) =
  gamma(constant(x)) + polyint(polyder(x)*gamma(restrict(x))*polygamma(0, restrict(x)))

_polygamma(k, x) = polygamma(k, constant(x)) + polyint(polyder(x)*polygamma(k + 1, restrict(x)))

# First two defs are to disambiguate with definitions in base
polygamma(k::Integer, x::AbstractSeries) = _polygamma(k, x)
polygamma{T<:Number}(a::AbstractArray{T}, x::AbstractSeries) =
  reshape([_polygamma(a[i], x) for i=1:length(a)], size(a))
polygamma(k, x::AbstractSeries) = _polygamma(k, x)

# TODO, what about the jump points
floor(x::AbstractSeries) = floor(constant(x)) + polyint(polyder(x)*0)
ceil(x::AbstractSeries) = ceil(constant(x)) + polyint(polyder(x)*0)
round(x::AbstractSeries) = round(constant(x)) + polyint(polyder(x)*0)
sign(x::AbstractSeries) = sign(constant(x)) + polyint(polyder(x)*0)
abs(x::AbstractSeries) = abs(constant(x)) + polyint(polyder(x)*sign(restrict(x)))

<(x::AbstractSeries, c::Real) = constant(x) < c
<(c::Real, x::AbstractSeries) = c < constant(x)
<(x::AbstractSeries, y::AbstractSeries) = constant(x) < constant(y)

generate(7)

export series, restrict, constant, polyint, polyval, polyder

end
