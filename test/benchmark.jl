using PowerSeries

function series_loop()
  accum = 0.0
  for i = 1:10000
    accum += sin(Series(1.0*i, 1.0)).ep
  end
end

function dual_loop()
  accum = 0.0
  for i = 1:10000
    accum += sin(Dual(1.0*i, 1.0)).ep
  end
end

function bare_loop()
  accum = 0.0
  for i = 1:10000
    accum += cos(1.0*i)
  end
end

# Warm up JIT
@assert series_loop() == bare_loop()
@assert dual_loop() == bare_loop()

tseries = @elapsed series_loop()
tdual = @elapsed dual_loop()
tbare = @elapsed bare_loop()

println("Series time")
println("Series", ", ", "Bare", ", ", "Ratio")
println(tseries, ", ", tbare, ", ", tseries/tbare)

println("\n")
println("Dual time")
println("Dual", ", ", "Bare", ", ", "Ratio")
println(tdual, ", ", tbare, ", ", tdual/tbare)

println("\n")
println("Series code")
code_native(series_loop, ())

println("\n")
println("Dual code")
code_native(dual_loop, ())

println("\n")
println("Bare code")
code_native(bare_loop, ())