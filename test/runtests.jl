using PowerSeries

# Represent the series 1.0 + 1.0*x - 2.0*x^2
a = series(1.0, 1.0, -2.0)

# Represent the series 1.0 + 0.0*x + 1.0*x^2
b = series(1.0, 0.0, 1.0)

# Series add linearly
@assert a+b == series(2.0, 1.0, -1.0)

# The output of series operations is truncated to match the input size.
# Represents (1+x-2x^2)(1 + x^2) = 1+x-x^2+o(x^3)
@assert a*b == series(1.0,1.0,-1.0)

# Functions with known derivatives can easily be overloaded to operate on
# power series.
# You can generate the taylor series of a function about a point x up to
# e.g. 6th order by computing f(Series(x, 1.0, 0.0, 0.0, 0.0, 0.0))
@assert sin(series(0.0, 1.0, 0.0, 0.0, 0.0, 0.0)) == series(0.0, 1.0, 0.0, -1.0/6.0, 0.0, 1.0/120)
@assert log(series(1.0, 1.0, 0.0, 0.0, 0.0, 0.0)) == series(0.0, 1.0, -1.0/2, 1.0/3, -1.0/4, 1.0/5)