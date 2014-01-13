using PowerSeries

@assert sin(series(0.0, 1.0, 0.0, 0.0)) == series(0.0, 1.0, 0.0, -1.0/6)
@assert log(series(1.0, 1.0, 0.0, 0.0)) == series(0.0, 1.0, -1.0/2, 1.0/3)