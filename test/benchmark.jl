using PowerSeries

function f()
  accum = 0.0
  for i = 1:10000
    accum += sin(Series(1.0*i, 1.0)).ep
  end
end

function g()
  accum = 0.0
  for i = 1:10000
    accum += cos(1.0*i)
  end
end

# Warm up JIT
f()
g()

t = @elapsed f()
tbase = @elapsed g()

print(t, ", ", tbase, ", ", t/tbase)