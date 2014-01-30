using PowerSeries


function series_loop()
  accum = 0.0
  for i = 1:10000
    accum += polyder(sin(series(1.0*i, 1.0)))
  end
  return accum
end

function bare_loop()
  accum = 0.0
  for i = 1:10000
    accum += cos(1.0*i)
  end
  return accum
end

# Warm up JIT
@assert series_loop() == bare_loop()

tseries = @elapsed series_loop()
tbare = @elapsed bare_loop()

println("Series time")
println("Series", ", ", "Bare", ", ", "Ratio")
println(tseries, ", ", tbare, ", ", tseries/tbare)
