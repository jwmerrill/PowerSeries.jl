typ(n) = symbol(string("Series", n))
elt(n) = symbol(string("c", n))
mem(s,n) = Expr(:., s, Expr(:quote, elt(n)))

function generate_type(n::Integer)
  @assert n > 0
  LastType = typ(n-1)
  Typ = typ(n)

  defn = :(immutable $Typ{T} <: AbstractSeries{T, $n} c0::T end)

  for i = 1:n
    push!(defn.args[3].args, :($(elt(i))::T))
  end

  eval(defn)

  if (n > 1)
    @eval restrict(x::$Typ) = $LastType($([mem(:x, i) for i = 0:n-1]...))
    @eval polyder(x::$Typ) = $LastType($([:($i*$(mem(:x, i))) for i = 1:n]...))
    @eval polyint{T}(x::$LastType{T}) = $Typ(zero(T),
      $([:($(mem(:x, i-1))/$i) for i = 1:n]...)
    )
  else
    @eval restrict(x::$Typ) = x.c0
    @eval polyder(x::$Typ) = x.c1
    @eval polyint(x::Real) = Series1(zero(typeof(x)), x)
  end

  #Horner's method
  polyval_body = mem(:x, n)
  for m = n-1:-1:0
    polyval_body = :($(mem(:x, m)) + eps*$polyval_body)
  end
  @eval polyval(x::$Typ, eps::Real) = polyval_body

  @eval constant(x::$Typ) = x.c0

  @eval +(x::$Typ, y::$Typ) = $Typ(
    $([:($(mem(:x, i)) + $(mem(:y, i))) for i = 0:n]...)
  )
  @eval +{T<:Real, S}(c::T, x::$Typ{S}) = $Typ{promote_type(T, S)}(
    c + $(mem(:x, 0)),
    $([mem(:x, i) for i = 1:n]...)
  )
  @eval +{T<:Real, S}(x::$Typ{S}, c::T) = $Typ{promote_type(T, S)}(
    $(mem(:x, 0)) + c,
    $([mem(:x, i) for i = 1:n]...)
  )

  @eval -(x::$Typ) = $Typ($([:(-$(mem(:x, i))) for i = 0:n]...))

  @eval -(x::$Typ, y::$Typ) = $Typ(
    $([:($(mem(:x, i)) - $(mem(:y, i))) for i = 0:n]...)
  )
  @eval -{T<:Real, S}(c::T, x::$Typ{S}) = $Typ{promote_type(T, S)}(
    c - $(mem(:x, 0)),
    $([:(-$(mem(:x, i))) for i = 1:n]...)
  )
  @eval -{T<:Real, S}(x::$Typ{S}, c::T) = $Typ{promote_type(T, S)}(
    $(mem(:x, 0)) - c,
    $([:(-$(mem(:x, i))) for i = 1:n]...)
  )

  times_term(i) = :(+($([:($(mem(:x, j))*$(mem(:y, i-j))) for j = 0:i]...)))
  @eval *(x::$Typ, y::$Typ) = $Typ($([times_term(i) for i = 0:n]...))
  @eval *(c::Real, x::$Typ) = $Typ($([:(c*$(mem(:x, i))) for i = 0:n]...))
  @eval *(x::$Typ, c::Real) = $Typ($([:($(mem(:x, i))*c) for i = 0:n]...))

  @eval /(c::Real, x::$Typ) = $Typ($([:(c/$(mem(:x, i))) for i = 0:n]...))
  @eval /(x::$Typ, c::Real) = $Typ($([:($(mem(:x, i))/c) for i = 0:n]...))
  
  # TODO probably a better way than eval to get a type from a string
  return eval(Typ)
end