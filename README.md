PowerSeries.jl
==============

Truncated Power Series for Julia

[![Build Status](https://travis-ci.org/jwmerrill/PowerSeries.jl.png?branch=master)](https://travis-ci.org/jwmerrill/PowerSeries.jl)

PowerSeries.jl defines Series types that represent truncated power series by their coefficients. You can do arithmetic on Series and apply functions to series just as you would Real or Complex numbers. Here's an example session:

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

# Extract the constant term of a series
julia> constant(a)
1.0

# Functions with known derivatives can easily be overloaded to operate on
# power series.
# You can generate the taylor series of a function about a point x up to
# e.g. 5th order by computing f(Series(x, 1.0, 0.0, 0.0, 0.0, 0.0))
julia> x = series(0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
Series5{Float64}(0.0,1.0,0.0,0.0,0.0,0.0)

julia> sin(x)
Series5{Float64}(0.0,1.0,0.0,-0.16666666666666666,0.0,0.008333333333333333)

julia> log(1 + x)
Series5{Float64}(0.0,1.0,-0.5,0.3333333333333333,-0.25,0.2)

julia> 1/(1 - x)
Series5{Float64}(1.0,1.0,1.0,1.0,1.0,1.0)

# These are numerically equal to their series definitions
julia> series(0.0, 1.0, 0.0, -1.0/6.0, 0.0, 1.0/120)
Series5{Float64}(0.0,1.0,0.0,-0.16666666666666666,0.0,0.008333333333333333)

julia> series(0.0, 1.0, -1.0/2, 1.0/3, -1.0/4, 1.0/5)
Series5{Float64}(0.0,1.0,-0.5,0.3333333333333333,-0.25,0.2)

# Take the derivative of a series
julia> polyder(a)
Series1{Float64}(1.0,-4.0)

# Integrate a series term by term. Note that by convention, the constant term is 0.
julia> polyint(a)
Series3{Float64}(0.0,1.0,0.5,-0.6666666666666666)

julia> @assert polyder(polyint(a)) == a

# Truncate a series to a series 1 order lower
julia> restrict(a)
Series1{Float64}(1.0,1.0)

# Restricting a first order series returns a real number
julia> restrict(restrict(a))
1.0

# polyint, polydir, and restrict are the only operations that change the order of
# a series. Arithmetic on series of different orders is disallowed because
# relevant terms in the lower order series may have been dropped at intermediate
# steps.
julia> series(1.0, 1.0) + series(1.0, 2.0, 3.0)
ERROR: no promotion exists for Series1{Float64} and Series2{Float64}
 in + at promotion.jl:158

# Truncated power series offer one of the best ways to take multiple derivatives
# of generic mathematical functions.
julia> f(x) = exp(-x^2)
f (generic function with 1 method)

julia> f2(x) = polyder(polyder(f(series(x, 1, 0))))
f2 (generic function with 1 method)

julia> f2(2.0)
0.25641894444227853

# Compare to the symbolic second derivative
julia> let x = 2.0; -2exp(-x^2)+4x^2*exp(-x^2); end
0.25641894444227853

# PowerSeries comes with types defined for series up to order 7. By default,
# trying to construct a higher order series is a type error.
julia> series(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
ERROR: no method series(Int64, Int64, Int64, Int64, Int64, Int64, Int64, Int64, Int64, Int64)

# If you want to work with higher order series, you can generate types up
# to a given order with PowerSeries.generate(order)
julia> PowerSeries.generate(9)
julia> series(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
Series9{Int64}(0,1,2,3,4,5,6,7,8,9)
```

For taking first derivatives of code, see also [DualNumbers.jl](https://github.com/scidom/DualNumbers.jl), and for taking symbolic derivatives, see the `differentiate` method of [Calculus.jl](https://github.com/johnmyleswhite/Calculus.jl).

Truncated series have performance advantages over symbolic derivatives for either deeply nested functions or high order derivatives.

###Theory of operation
Computations of functions of a power series are based on the fundamental theorem of calculus:

![equation-1](http://latex.codecogs.com/png.latex?f%28x%20+%20%5Cepsilon%29%20%3D%20f%28x%29%20+%20%5Cint_x%5E%7Bx%20+%20%5Cepsilon%7D%20dx%20f%27%28x%29)

Using this relation, it's easy to derive a composition rule for functions that can be applied directly to power series.

![equation-2](http://latex.codecogs.com/png.latex?f%28g%28x%20+%20%5Cepsilon%29%29%20%3D%20f%28g%28x%29%29%20+%20%5Cint_x%5E%7Bx%20+%20%5Cepsilon%7D%20dx%20f%27%28g%28x%29%29%20g%27%28x%29)

This is essentially an extension of the chain rule from infinitesimal calculus to finite step sizes.

Once differentiation and definite integration are defined on series, this relation allows a simple definition of functions of series. For example, the sine and cosine of series are mutually-recursively defined as

```julia
sin(x::AbstractSeries) = sin(constant(x)) + polyint(polyder(x)*cos(restrict(x)))
cos(x::AbstractSeries) = cos(constant(x)) - polyint(polyder(x)*sin(restrict(x)))
```

The general pattern is

```julia
f(x::AbstractSeries) = f(constant(x)) + polyint(polyder(x)*f'(restrict(x)))
```

where `f'` should be replaced by a known derivative function.
