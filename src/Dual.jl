import Base: sin

immutable Dual{T<:Number} <: Number
  re::T
  ep::T
end

sin(p::Dual) = Dual(sin(p.re), p.ep*cos(p.re))

export Dual