module PowerSeries

import Base: sin, cos, exp, log

abstract AbstractSeries{T<:Real, N} <: Number

include("generate_type.jl")

series_types = DataType[]

# Generate series types up to 3rd order
push!(series_types, generate_type(1))
push!(series_types, generate_type(2))
push!(series_types, generate_type(3))

function series(n...)
  l = length(n)

  # Generate new types on demand
  while l > length(series_types)
    push!(series_types, generate_type(length(series_types) + 1))
  end

  return series_types[l - 1](promote(n...)...)
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