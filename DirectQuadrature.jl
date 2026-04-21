module DirectQuadrature

using SpecialFunctions


#TODO: Add comment and citation to InverseLaplace.jl
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

#=
function DQ_Array(input_func_array, t, R, N; γ = 0)
    Dt = typeof(t)
    Dt = Dt <: Complex ? Dt : Complex{Dt}
    if Dt <: Int
      @error "Type conversions fail for integers. Consider using a t with type BigFloat."
    end


  sVals = γ .+  1im * LinRange(-R, R, N)
  Δs = sVals[2] - sVals[1]
  return Δs * (1/( 2*convert(Dt, π)*1im ) )*sum(input_func_array.(sVals) .* exp.(sVals * t))


end
=#


function  DQ_Array(input_func_array, t, R, N; γ = 0)
  DQ(input_func_array, t, R, N; γ = γ)
end




end
