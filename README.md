PowerSeries.jl
==============

Truncated Power Series for Julia

PowerSeries.jl exports a Series type that represents a truncated power series by its coefficients. You can do arithmetic on Series and apply functions to series just as you would Real or Complex numbers. Here's an example session:

```julia
julia> using PowerSeries

# Represent the series 1.0 + 1.0*x - 2.0*x^2
julia> a = Series(1.0, 1.0, -2.0)
Series{Float64}(1.0,1.0,-2.0)

# Represent the series 1.0 + 0.0*x + 1.0*x^2
julia> b = Series(1.0, 0.0, 1.0)
Series{Float64}(1.0,0.0,1.0)

# Series add linearly
julia> a+b
Series{Float64}(2.0,1.0,-1.0)

# The output of series operations is truncated to match the input size.
# Represents (1+x-2x^2)(1 + x^2) = 1+x-x^2+o(x^3)
julia> a*b
Series{Float64}(1.0,1.0,-1.0)

# Functions with known derivatives can easily be overloaded to operate on
# power series.
# You can generate the taylor series of a function about a point x up to
# e.g. 6th order by computing f(Series(x, 1.0, 0.0, 0.0, 0.0, 0.0))
julia> sin(Series(0.0, 1.0, 0.0, 0.0, 0.0, 0.0))
Series{Float64}(0.0,1.0,0.0,-0.16666666666666666,0.0,0.008333333333333333)
julia> log(Series(1.0, 1.0, 0.0, 0.0, 0.0, 0.0))
Series{Float64}(0.0,1.0,-0.5,0.3333333333333333,-0.25,0.2)

# These are numerically equal to their series definitions
julia> Series(0.0, 1.0, 0.0, -1.0/6.0, 0.0, 1.0/120)
Series{Float64}(0.0,1.0,0.0,-0.16666666666666666,0.0,0.008333333333333333)
julia> Series(0.0, 1.0, -1.0/2, 1.0/3, -1.0/4, 1.0/5)
Series{Float64}(0.0,1.0,-0.5,0.3333333333333333,-0.25,0.2)
```

###Status
PowerSeries.jl is currently a proof of concept. The API may change before it is released as a Julia package. The internal representation of Series is _very_ likely to change to improve performance.

###Theory of operation
Computations of functions of a power series are based on the fundamental theorem of calculus:

```
f(x + \epsilon) = f(x) + \int_0^\epsilon dx f'(x + \epsilon)
```

Using this relation, it's easy to derive a composition rule for functions that can be applied directly to power series.

```
f(g(x + \epsilon)) = f(g(x)) + \int_0^\epsilon dx f'(g(x + \epsilon)) g'(x + \epsilon)
```

This is essentially an extension of the chain rule from infinitesimal calculus to finite step sizes.

Once differentiation and definite integration are defined on series, this relation allows a simple definition of functions of series. For example, the sin of a series is defined as

```julia
sin(p::Series) = sin(p.re) + pint(diff(p)*cos(restrict(p)))
```

where `p.re` is the "real" or "constant" part of the series `p`, `pint` is the definite integral, `diff` is the series derivative, and `restrict` truncates a series of order `o(x^n)` to `o(x^{n-1})`.

###Internal representation
Series are currently represented as linked lists with a base type number as the head, and a series over the same type as the tail.

```julia
immutable Series{T<:Number} <: Number
  re::T
  ep::Union(T, Series{T})
end
```

This allows simple recursive definitions of all operations, but appears to carry a significant performance overhead. I'm currently evaluating switching to a tower of fixed length types, or a type backed by an array to improve performance.