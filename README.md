PowerSeries.jl
==============

Truncated Power Series for Julia

PowerSeries.jl exports Series types that represent truncated power series by their coefficients. You can do arithmetic on Series and apply functions to series just as you would Real or Complex numbers. Here's an example session:

```julia
julia> using PowerSeries

# Represent the series 1.0 + 1.0*x - 2.0*x^2
julia> a = series(1.0, 1.0, -2.0)
Series2{Float64}(1.0,1.0,-2.0)

# Represent the series 1.0 + 0.0*x + 1.0*x^2
julia> b = series(1.0, 0.0, 1.0)
Series2{Float64}(1.0,0.0,1.0)

# Series add linearly
julia> a+b
Series2{Float64}(2.0,1.0,-1.0)

# The output of series operations is truncated to match the input size.
# Represents (1+x-2x^2)(1 + x^2) = 1+x-x^2+o(x^3)
julia> a*b
Series2{Float64}(1.0,1.0,-1.0)

# Functions with known derivatives can easily be overloaded to operate on
# power series.
# You can generate the taylor series of a function about a point x up to
# e.g. 5th order by computing f(Series(x, 1.0, 0.0, 0.0, 0.0, 0.0))
julia> sin(series(0.0, 1.0, 0.0, 0.0, 0.0, 0.0))
Series5{Float64}(0.0,1.0,0.0,-0.16666666666666666,0.0,0.008333333333333333)
julia> log(series(1.0, 1.0, 0.0, 0.0, 0.0, 0.0))
Series5{Float64}(0.0,1.0,-0.5,0.3333333333333333,-0.25,0.2)

# These are numerically equal to their series definitions
julia> series(0.0, 1.0, 0.0, -1.0/6.0, 0.0, 1.0/120)
Series5{Float64}(0.0,1.0,0.0,-0.16666666666666666,0.0,0.008333333333333333)
julia> series(0.0, 1.0, -1.0/2, 1.0/3, -1.0/4, 1.0/5)
Series5{Float64}(0.0,1.0,-0.5,0.3333333333333333,-0.25,0.2)
```

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
sin(x::AbstractSeries) = sin(constant(x)) + polyint(polyder(x)*cos(restrict(x)))
```

where `constant(x)` is the constant term in the series `x`, `polyint` is the definite integral, `polyder` is the series derivative, and `restrict` truncates a series of order `o(x^n)` to `o(x^{n-1})`.
