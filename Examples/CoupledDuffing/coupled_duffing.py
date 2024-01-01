from DQbee import *
import sympy as sp
import pandas as pd
import sys
import time
sys.path.append("../../")


# Defining the coupling and dumping parameters
def A_entry(i, j):
    if i == j:
        return 1
    if i - j in [-1, 1]:
        return sp.Rational(1, 3)
    return 0


def compute_coupled_duffing(N, methods):
    # Creating variables for positions and velocities
    xs = sp.symbols(" ".join([f"x{i}" for i in range(N)]))
    dxs = sp.symbols(" ".join([f"dx{i}" for i in range(N)]))

    A = sp.Matrix([[A_entry(i, j) for j in range(N)] for i in range(N)])
    delta = 2

    # Defining the equations
    eqs = None
    if N == 1:
        eqs = [sp.Eq(xs, dxs)]
        linear = A[0] * xs
        eqs.append(sp.Eq(dxs, linear - linear**3 - delta * dxs))
    else:
        eqs = [sp.Eq(xs[i], dxs[i]) for i in range(N)]
        for i in range(N):
            linear = (A * sp.Matrix(xs))[i]
            eqs.append(sp.Eq(dxs[i], linear - linear**3 - delta * dxs[i]))
    coupled_duffing = EquationSystem(eqs)

    # Computing inner quadratization only
    start = time.time()
    inner_quadratization = optimal_inner_quadratization(
        coupled_duffing, display=False)
    num_new_variables = len(inner_quadratization[3])
    end = time.time()
    inner_quadratic_time = end - start
    print(f"Inner quadratization computed in {end - start} seconds")

    # Computing the dissipative equilibria
    equilibria = []
    for rhs in range(2**N):
        rhs_vect = []
        for i in range(N):
            if rhs % 2 == 1:
                rhs_vect.append(1)
            else:
                rhs_vect.append(-1)
            rhs = rhs // 2
        equilibria.append(
            list(A.solve(sp.Matrix(rhs_vect))) + [0 for _ in range(N)])

    # Computing the dissipative quadratization
    for method in methods:  # ['numpy', 'Routh-Hurwitz']:
        start = time.time()
        dissipative_multi_system = dquadratization_multi_equilibrium(
            coupled_duffing,
            equilibria,
            method=method,
            display=False
        )
        end = time.time()
        dissipative_time = (end - start, method)
        lambda_value = dissipative_multi_system[0]
        print(
            f"Diffipative quadratization computed in {end - start} seconds ({method} method)")

    return inner_quadratic_time, dissipative_time, num_new_variables, lambda_value


def main():
    methods = ['numpy']
    results = []

    for N in range(1, 9):
        inner_quadratic_time, dissipative_time, num_new_variables, lambda_value = compute_coupled_duffing(
            N, methods)
        results.append([N, num_new_variables, inner_quadratic_time,
                        dissipative_time[0], dissipative_time[1], lambda_value])

    df = pd.DataFrame(results, columns=[
        "N", "Number of new vars", "Time (inner-quadratic)", "Time (dissipative)", "Method", "Lambda"])
    df.to_excel("results_numpy.xlsx", index=False)

# def main():
#     methods = ['Routh-Hurwitz']
#     results = []

#     for N in range(1, 5):
#         inner_quadratic_time, dissipative_time, num_new_variables, lambda_value = compute_coupled_duffing(
#             N, methods)
#         results.append([N, num_new_variables, inner_quadratic_time,
#                         dissipative_time[0], dissipative_time[1], lambda_value])

#     df = pd.DataFrame(results, columns=[
#         "N", "Number of new vars", "Time (inner-quadratic)", "Time (dissipative)", "Method", "Lambda"])
#     df.to_excel("results_routh.xlsx", index=False)


main()
