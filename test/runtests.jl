using PowerSeries

@assert sin(Series(0.0, 1.0, 0.0, 0.0, 0.0, 0.0)) == Series(0.0, 1.0, 0.0, -1.0/6, 0.0, 1.0/120)
@assert log(Series(1.0, 1.0, 0.0, 0.0, 0.0, 0.0)) == Series(0.0, 1.0, -1.0/2, 1.0/3, -1.0/4, 1.0/5)