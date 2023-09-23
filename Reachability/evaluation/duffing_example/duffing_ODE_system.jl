# =================
# Dependencies
# =================

using ReachabilityAnalysis, CarlemanLinearization
using Plots, LaTeXStrings, LinearAlgebra, SparseArrays
using Plots.PlotMeasures
using LazySets: center
using CarlemanLinearization: _error_bound_specabs_Re_λ₁, _error_bound_specabs_R

include("../utils.jl")

# ====================================================================
# Model definition - The example of the Duffing oscillator
# ====================================================================

function duffing_system_carlin()
    F1 = zeros(3, 3)
    F1[1, 2] = 1
    F1[2, 1] = -1
    F1[2, 2] = -1
    F1[3, 3] = -1

    F2 = zeros(3, 9) # [x, x⊗x]
    F2[2, 3] = 1
    F2[3, 1] = 1
    F2[3, 2] = 2

    return F1, F2
end

# ====================================================================
# Solution method - with carleman linearization
# ====================================================================

function _solve_system_carlin(; N=4, T=30.0, δ=0.1, radius0=0, bloat=false, resets=nothing, x0=0.1, x1=0.1)
    x0c = [0.1, 0.1, 0.01]

    F1, F2 = duffing_system_carlin()
    Re_λ₁ = _error_bound_specabs_Re_λ₁(F1; check=false)
    R = _error_bound_specabs_R(x0c, F2, Re_λ₁)
    println("R = ", R)
    println("Re_λ₁ = ", Re_λ₁)

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

# ====================================================================
# Plot the graph
# ====================================================================

# Ploting the results
Tmax = 20.0
rr0 = 0.0

sol_N4_no_error = _solve_system_carlin(N=6, T=Tmax, δ=0.02, radius0=rr0, bloat=false)
sol_N4_error = _solve_system_carlin(N=6, T=Tmax, δ=0.02, radius0=0, bloat=true)

# figure with NO error bounds, plot for x
function figure_System_NoError()

    fig = plot(legend=:topright, 
                xlab = L"\textrm{Time t}", 
                ylab = L"\textrm{x(t)} ",
                legendfontsize=12,
                tickfont=font(14, "Times"),  
                guidefontsize=14,  
                xtick = 0:5:20,
                xguidefont=font(14, "Times"),  
                yguidefont=font(14, "Times"),  
                bottom_margin=5mm,
                left_margin=5mm,
                right_margin=5mm,
                top_margin=5mm,
                size=(800, 600),
                grid = true) 

    plot!(fig, sol_N4_no_error,  
          vars=(0, 1), 
          color=:aquamarine, 
          lc=:aquamarine, 
          linewidth=3,
          linestyle=:dash, 
          label=L"x''=-x+x^{3}-x'")

    return fig

end


fig = figure_System_NoError()
savefig(fig, joinpath(TARGET_FOLDER, "figure_duffing_no_error.pdf"))
display(fig)

# figure with error bounds, plot for x
function figure_System_withError()

    fig = plot(legend=:topright, 
                xlab = L"\textrm{Time t}", 
                ylab = L"\textrm{x(t)} ",
                legendfontsize=12,
                tickfont=font(14, "Times"),  
                guidefontsize=14,  
                xguidefont=font(14, "Times"),  
                yguidefont=font(14, "Times"),  
                bottom_margin=5mm,
                left_margin=5mm,
                right_margin=5mm,
                top_margin=5mm,
                size=(800, 600),
                grid = true) 

    plot!(fig, sol_N4_error,  
            vars=(0, 1), 
            color=:orange, 
            lc=:orange, 
            linewidth=3,
            linestyle=:dash, 
            label=L"\text{error bounds, } N=6")

    plot!(fig, sol_N4_no_error,  
          vars=(0, 1), 
          color=:aquamarine, 
          lc=:aquamarine, 
          linewidth=3,
          linestyle=:dash, 
          label=L"x''=-x+x^{3}-x'")

    return fig

end

fig = figure_System_withError()
savefig(fig, joinpath(TARGET_FOLDER, "figure_duffing_error.pdf"))
display(fig)