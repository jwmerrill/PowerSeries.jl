using Base.Test
using PowerSeries

# Represent the series 1.0 + 1.0*x - 2.0*x^2
a = series(1.0, 1.0, -2.0)

# Represent the series 1.0 + 0.0*x + 1.0*x^2
b = series(1.0, 0.0, 1.0)

# Series add linearly
@test a+b == series(2.0, 1.0, -1.0)

# The output of series operations is truncated to match the input size.
# Represents (1+x-2x^2)(1 + x^2) = 1+x-x^2+o(x^3)
@test a*b == series(1.0,1.0,-1.0)

# Functions with known derivatives can easily be overloaded to operate on
# power series.
# You can generate the taylor series of a function about a point x up to
# e.g. 6th order by computing f(Series(x, 1.0, 0.0, 0.0, 0.0, 0.0))
x = series(0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
@test sin(x) == series(0.0, 1.0, 0.0, -1.0/6.0, 0.0, 1.0/120)
@test log(1 + x) == series(0.0, 1.0, -1.0/2, 1.0/3, -1.0/4, 1.0/5)
@test 1/(1 - x) == series(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

@test polyval(a, 0.1) == 1.08
@test polyval(a, -0.1) == 0.88
@test polyder(a) == series(1.0, -4.0)
@test polyint(a) == series(0.0, 1.0, 1.0/2, -2.0/3)
@test restrict(a) == series(1.0, 1.0)

# Power series are one of the best ways to take higher order derivatives
# of generic functions.
#
# This example tests promotion rules in the constructor because of the
# 1 and 0, instead of 1.0 and 0.0
f(x) = exp(-x^2)
f2(x) = polyder(polyder(f(series(x, 1, 0))))
@test f2(2.0) == let x = 2.0; -2exp(-x^2)+4x^2*exp(-x^2); end

# macro-expansion-time conditional.
# This is probably not very good or general, but gets around the fact that
# the signture for @test_throws changed between versions 0.2 and 0.3
macro when(cond, ex)
  eval(cond) ? esc(ex) : nothing
end

# PowerSeries comes with types defined for series up to order 7. By default,
# trying to construct a higher order series is a type error.
@when(
  VERSION >= VersionNumber(0,2,typemax(Int64)),
  @test_throws MethodError series(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
)

# If you want to work with higher order series, you can generate types up
# to a given order with PowerSeries.generate(order)
PowerSeries.generate(9)
@test series(0, 1, 2, 3, 4, 5, 6, 7, 8, 9) != false

@test series(1, 2) - 1 == series(0, 2)

# Regression test for https://github.com/jwmerrill/PowerSeries.jl/issues/6
let x=series(0.0, 1.0, 0.0, 0.0)
  @test x*x == x^2
end

@test series(1.0, 2.0, 3.0) < 2.0
@test 2.0 > series(1.0, 2.0, 3.0)
@test 0.5 < series(1.0, 2.0, 3.0)
@test series(1.0, 2.0, 3.0) > 0.5
@test !(series(1.0, 2.0, 3.0) < 1.0)
@test !(series(1.0, 2.0, 3.0) > 1.0)
@test series(1.0, 1.0) < series(2.0, 0.0)
@test series(2.0, 0.0) > series(1.0, 1.0)
@test !(series(1.0, 2.0) < series(1.0, 0.0))
@test !(series(1.0, 0.0) > series(1.0, 2.0))

# Only 1 argument functions
fns = [
  sqrt,
  exp,
  log,
  sin,
  cos,
  tan,
  asin,
  acos,
  atan,
  sinh,
  cosh,
  tanh,
  asinh,
  acosh,
  atanh,
  csc,
  sec,
  cot,
  acsc,
  asec,
  acot,
  csch,
  sech,
  coth,
  acsch,
  asech,
  acoth,
  gamma,
  floor,
  ceil,
  round,
  sign,
  abs
]

for val in [-2.3, -1.2, -0.1, 0.7, 3.6]
  for fn in fns
    if (val < 1 && fn === acosh) continue end
    if (val < 0 && fn === sqrt) continue end
    if (val < 0 && fn === log) continue end
    if (abs(val) > 1 && fn === asin) continue end
    if (abs(val) > 1 && fn === acos) continue end
    if (abs(val) > 1 && fn === atanh) continue end
    if ((val < 0 || val > 1) && fn === asech) continue end
    if (abs(val) < 1 && fn === acsc) continue end
    if (abs(val) < 1 && fn === asec) continue end
    if (abs(val) < 1 && fn === acoth) continue end
    @test_approx_eq_eps constant(polyder(fn(series(val, 1.0)))) 5000*(fn(val + .0001) - fn(val- .0001)) .01
  end
end
