# =================
# Dependencies
# =================

using ReachabilityAnalysis, CarlemanLinearization
using Plots, LaTeXStrings, LinearAlgebra, SparseArrays
using Plots.PlotMeasures
using LazySets: center
using CarlemanLinearization: _error_bound_specabs_R

include("../utils.jl")

# =================
# Model definition (The second example)
# =================

function system_carlin(a, k, m)

  F1 = zeros(3, 3)
  F1[1, 2] = 1
  F1[2, 1] = k
  F1[3, 2] = -m

  F2 = zeros(3, 9) # [x, x⊗x]
  F2[2, 3] = a
  F2[3, 2] = 2

  return F1, F2
end

# See Definition (2.2) in [2]. These bounds use the spectral norm (p = 2)
function _error_bound_specabs_R_(x₀, F₁, F₂; check=true)
  nx₀ = norm(x₀, 2)
  nF₂ = opnorm(F₂, 2)

  # compute eigenvalues and sort them by increasing real part
  λ = eigvals(F₁, sortby=real)
  println("λ = ", λ)
  λ₁ = last(λ)
  Re_λ₁ = real(λ₁)
  print("The real part of the largest eigenvalue is ", Re_λ₁)
  if check
      @assert Re_λ₁ <= 0 "expected Re(λ₁) ≤ 0, got $Re_λ₁"
  end
  R = nx₀ * nF₂ / abs(Re_λ₁)
  return (R, Re_λ₁)
end

function _solve_system_carlin(; N=4, T=30.0, δ=0.1, radius0=0, bloat=false, resets=nothing, x0, x1, a, m, k)
  x0c = [x0, x1, x0^2 - m * x0]

  F1, F2 = system_carlin(a, k, m)
  R, Re_lambda1 = _error_bound_specabs_R_(x0c, F1, F2; check=true)

  n = 3
  dirs = _template(n=n, N=N)
  alg = LGG09(δ=δ, template=dirs)

  if radius0 == 0
    X0 = convert(Hyperrectangle, Singleton(x0c))
  else
    X0 = Hyperrectangle(x0c, radius0)
  end

  if isnothing(resets)
    @time sol = _solve_CARLIN(X0, F1, F2; alg=alg, N=N, T=T, bloat=bloat)
  else
    @time sol = _solve_CARLIN_resets(X0, F1, F2; resets=resets, alg=alg, N=N, T=T, bloat=bloat)
  end

  return sol
end

# @taylorize function system_equation(dx, x, a, k, m)
#   x1, x2, x3 = x # y = x2

#   dx[1] = x2
#   dx[2] = k * x1 + a * x1 * x3
#   dx[3] = 2 * x1 * x2 - m * x2

# end

# Ploting the results
Tmax = 10.0
rr0 = 0.0

# figure with error bounds with comparison
solN4_a1_bloat = _solve_system_carlin(N=4, T=Tmax, δ=0.01, radius0=rr0, bloat=true, x0 = 0.1, x1 = 0.1, a = 1.0, m = 1.0, k = - 3.0)