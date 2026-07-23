module DirectQuadrature

using SpecialFunctions


"""
    DQ(func, t, R, N; γ = 0)

Return the inverse Laplace transform of `func` evaluated at a time `t`, by integrating between `γ - iR` and `γ + iR` using a Riemann sum with `N` quadrature points.
Note that `γ` should be greater than the real component of all poles of `func`.
"""
function DQ(func, t, R, N; γ = 0)

  Dt = typeof(t)
  Dt = Dt <: Complex ? Dt : Complex{Dt}
  if Dt <: Int
      @error "Type conversions fail for integers. Consider using a t with type BigFloat."
  end

  sVals = γ .+  1im * LinRange(-R, R, N)
  Δs = sVals[2] - sVals[1]
  return Δs * (1/( 2*convert(Dt, π)*1im ) )*sum(func.(sVals) .* exp.(sVals * t))
end


"""
    DQ_Array(input_func_array, t, R, N; γ = 0)

Return the inverse Laplace transform of `input_func_array` evaluated at a time `t`, by integrating between `γ - iR` and `γ + iR` using `N` quadrature points.
Note that `γ` should be greater than the real component of all poles of `input_func_array`.

This is an alias for `DQ` which handles arrays by default. It is provided for consistency.
"""
function  DQ_Array(input_func_array, t, R, N; γ = 0)
  DQ(input_func_array, t, R, N; γ = γ)
end




end
