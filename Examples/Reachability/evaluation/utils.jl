using ReachabilityAnalysis
using LazySets.Arrays

# functns for Kronecker powers and sums
using ReachabilityAnalysis: kron_pow, kron_pow_stack

using LaTeXStrings

__toL(x, digits) = L"%$(round(x, digits=digits))"
xticklatex(vec, digits) = (vec, __toL.(vec, Ref(digits)));

# functions for error bounds
using CarlemanLinearization
using CarlemanLinearization: kron_sum,
    error_bound_pseries,
    error_bound_specabs,
    convergence_radius_specabs,
    _error_bound_specabs_R,
    build_matrix

using ReachabilityAnalysis: ReachSolution

# TODO remove after ReachabilityAnalysis#551
LazySets.Approximations.box_approximation(X::Vector{IntervalArithmetic.Interval{Float64}}) = box_approximation(IntervalBox(X))

function _template(; n, N)
    dirs = Vector{SingleEntryVector{Float64}}()
    d = sum(n^i for i in 1:N)

    for i in 1:n
        x = SingleEntryVector(i, d, 1.0)
        push!(dirs, x)
    end
    for i in 1:n
        x = SingleEntryVector(i, d, -1.0)
        push!(dirs, x)
    end
    return CustomDirections(dirs)
end

function _project(sol, vars)
    πsol_1n = Flowpipe([ReachSet(set(project(R, vars)), tspan(R)) for R in sol])
end

# general method given a reachability algorithm for the linear system
function _solve_CARLIN_alg(X0, F1, F2, alg;
    resets=0, N, T, Δt0=interval(0),
    error_bound_func=error_bound_specabs)
    if resets == 0
        _solve_CARLIN(X0, F1, F2; alg=alg, N=N, T=T, Δt0=interval(0), error_bound_func=error_bound_func)
    else
        _solve_CARLIN_resets(X0, F1, F2; resets=resets, alg=alg, N=N, T=T, Δt0=interval(0), error_bound_func=error_bound_func)
    end
end

# method using the LGG09 reachability algorithm
function _solve_CARLIN_LGG09(X0, F1, F2;
    resets=0, N, T, δ, Δt0=interval(0),
    error_bound_func=error_bound_specabs)
    n = dim(X0)
    dirs = _template(n=n, N=N)
    alg = LGG09(δ=δ, template=dirs)
    return _solve_CARLIN_alg(X0, F1, F2, alg;
        resets=resets, N=N, T=T, Δt0=Δt0,
        error_bound_func=error_bound_func)
end

function _solve_CARLIN(X0, F1, F2; alg, N, T, Δt0=interval(0), A=nothing,
    error_bound_func=error_bound_specabs,
    bloat=true, kron_pow_algorithm="explicit")

    # lift initial states
    n = dim(X0)
    ŷ0 = kron_pow_stack(X0, N) |> box_approximation

    # solve continuous ODE
    if isnothing(A)
        A = build_matrix(F1, F2, N)
    end
    prob = @ivp(ŷ' = Aŷ, ŷ(0) ∈ ŷ0)
    sol = solve(prob, T=T, alg=alg, Δt0=Δt0)

    # projection onto the first n variables
    πsol_1n = _project(sol, 1:n)

    if !bloat
        return πsol_1n
    end

    # compute errors
    errfunc = error_bound_func(X0, Matrix(F1), Matrix(F2), N=N)

    # evaluate error bounds for each reach-set in the solution
    E = [errfunc.(tspan(R)) for R in sol]

    # if the interval is always > 0 then we can just take max(Ei)

    # symmetrize intervals
    E_rad = [symmetric_interval_hull(Interval(ei)) for ei in E]
    E_ball = [BallInf(zeros(n), max(high(ei)...)) for ei in E_rad]

    # sum the solution with the error
    fp_bloated = [ReachSet(set(Ri) ⊕ Ei, tspan(Ri)) for (Ri, Ei) in zip(πsol_1n, E_ball)] |> Flowpipe

    return ReachSolution(fp_bloated, alg)
end

function _compute_resets(resets::Int, T)
    return mince(0 .. T, resets + 1)
end

function _compute_resets(resets::Vector{Float64}, T)
    # assumes initial time is 0
    aux = vcat(0.0, resets, T)
    return [interval(aux[i], aux[i+1]) for i in 1:length(aux)-1]
end

function _solve_CARLIN_resets(X0, F1, F2; resets, alg, N, T, Δt0=interval(0),
    error_bound_func=error_bound_specabs,
    bloat=true)

    # build state matrix (remains unchanged upon resets)
    A = build_matrix(F1, F2, N)

    # time intervals to compute
    time_intervals = _compute_resets(resets, T)

    # compute until first chunk
    T1 = sup(first(time_intervals))
    sol_1 = _solve_CARLIN(X0, F1, F2; alg=alg, N=N, T=T1, Δt0=interval(0), A=A, error_bound_func=error_bound_func, bloat=bloat)

    # preallocate output flowpipe
    fp_1 = flowpipe(sol_1)
    out = Vector{typeof(fp_1)}()
    push!(out, fp_1)

    # approximate final reach-set
    Rlast = sol_1[end]
    X0 = box_approximation(set(Rlast))

    # compute remaining chunks
    for i in 2:length(time_intervals)
        T0 = T1
        Ti = sup(time_intervals[i])
        sol_i = _solve_CARLIN(X0, F1, F2; alg=alg, N=N, T=Ti - T0, Δt0=interval(T0), A=A, error_bound_func=error_bound_func, bloat=bloat)
        push!(out, flowpipe(sol_i))
        X0 = box_approximation(set(sol_i[end]))
        T1 = Ti
    end

    return ReachSolution(MixedFlowpipe(out), alg)
end
