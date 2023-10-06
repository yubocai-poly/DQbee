import sys
sys.path.append("../DQbee/")
# from Simulation import *
import sympy as sp
from DQbee.Combinations import *
from DQbee.EquationSystem import *
from DQbee.DQuadratization import *
import pytest
import math
import time


# ---------------- Examples setting ----------------
# Example Van der Pol Oscillator
x, y, z = sp.symbols("x, y z")
mu = sp.Symbol("mu")
w1 = sp.Symbol("w1")

system_vanderpol = [
    sp.Eq(y, z),
    sp.Eq(z, mu * (1 - y ** 2) * z - y)
]

eq_system_vanderpol = EquationSystem(system_vanderpol)
oiq_vanderpol, sub_oiq_vanderpol, monomial_to_quadratic_form, map_variables = optimal_inner_quadratization(
    eq_system_vanderpol)
Diq_vanderpol, F1, map_variables = optimal_dissipative_quadratization(
    eq_system_vanderpol, sub_oiq_vanderpol, monomial_to_quadratic_form, map_variables, [-3])


# ---------------- Testing for Quadratization Part ----------------
def test_input_system():
    assert eq_system_vanderpol.system == system_vanderpol


def test_optimal_inner_quadratization():
    test_oiq_vanderpol_system = [sp.Eq(y, z), sp.Eq(
        z, mu*z*(1 - y**2) - y), sp.Eq(y**2, 2*y*z)]
    test_sub_oiq_vanderpol_system = [sp.Eq(y, z), sp.Eq(
        z, -mu*w1*z + mu*z - y), sp.Eq(w1, 2*y*z)]
    assert oiq_vanderpol.system == test_oiq_vanderpol_system and sub_oiq_vanderpol.system == test_sub_oiq_vanderpol_system


def test_dissipative_inner_quadratization():
    test_F1 = sp.Matrix([[0, 1, 0], [-1, mu, 0], [0, 0, -3]])
    test_Diq_vanderpol_system = [sp.Eq(y, z), sp.Eq(
        z, -mu*w1*z + mu*z - y), sp.Eq(w1, -3*w1 + 3*y**2 + 2*y*z)]
    assert Diq_vanderpol.system == test_Diq_vanderpol_system and F1 == test_F1


# ---------------- Testing for Weakly nonlinearity ----------------
def test_compute_weakly_nonlinearity():
    mu = -1
    system_vanderpol_ = [
        sp.Eq(y, z),
        sp.Eq(z, mu * (1 - y ** 2) * z - y)
    ]
    eq_system_vanderpol_ = EquationSystem(system_vanderpol_)
    test_upper_bound = 2 * math.sqrt(13)
    test_F2 = sp.Matrix([[0, 0, 0, 0, 0, 0, 0, 0, 0], [
        0, 0, 0, 0, 0, 1, 0, 0, 0], [3, 2, 0, 0, 0, 0, 0, 0, 0]])
    oiq_vanderpol_, sub_oiq_vanderpol_, monomial_to_quadratic_form_, map_variables_ = optimal_inner_quadratization(
        eq_system_vanderpol_)
    Diq_vanderpol_, F1_, map_variables_ = optimal_dissipative_quadratization(
        eq_system_vanderpol, sub_oiq_vanderpol_, monomial_to_quadratic_form_, map_variables_, [-3])
    F2, upper_bound = compute_weakly_nonlinearity(Diq_vanderpol_)
    assert (test_upper_bound - upper_bound) < 1e-4 and F2 == test_F2


# ---------------- Testing for Simulation ----------------
# initial_state = {y: 0.1, z: 0}
# symbolic_args = {mu: 1}


# def test_simulation_original_system():
#     state_original = system_to_odeint(eq_system_vanderpol, t=[
#                                       0, 20, 1000], symbolic_args=symbolic_args, initial_state=initial_state)
#     assert state_original is not None


# def test_simulation_DQ_system():
#     state_dissipative = system_to_odeint(Diq_vanderpol, t=[
#                                          0, 20, 1000], symbolic_args=symbolic_args, initial_state=initial_state)
#     assert state_dissipative is not None


# ---------------- Testing for larger degree system  ----------------
x1, x2 = sp.symbols("x1, x2")
system_large_degree = [
    sp.Eq(x1, x2**7),
    sp.Eq(x2, x1**3),
]

eq_system_large_degree = EquationSystem(system_large_degree)


def test_computing_time():
    start_time = time.time()
    result = optimal_inner_quadratization(eq_system_large_degree)
    end_time = time.time()
    assert (end_time - start_time) < 45


# ---------------- Testing for Algorithm 2  ----------------
duffing_system_1 = [
    sp.Eq(x1, x2),
    sp.Eq(x2, - x1 + x1 ** 3 - x2)
]

duffing_eq_system_1 = EquationSystem(duffing_system_1)


def test_algorithm_2():
    duffing_eq_system_1_one_equilibrium = dquadratization_one_equilibrium(
        duffing_eq_system_1, {x1: 0, x2: 0}, display=True)
    assert duffing_eq_system_1_one_equilibrium is not None


# ---------------- Testing for Algorithm 3  ----------------
# ---------------- Testing for larger dimension cases ----------------
x = sp.Symbol('x')
num_list = list(range(1, 8))

rhs = - x
for i in num_list:
    rhs = rhs * (x - i)

system = [sp.Eq(x, rhs)]
multistable_eq_system = EquationSystem(system)

rhs_differential = sp.diff(rhs, x, 1).expand()
dissipative_equilibrium_list = []
for i in range(1, 8):
    if rhs_differential.subs(x, i) < 0:
        dissipative_equilibrium_list.append([i])


def test_algorithm_3_large_numpy():
    dissipative_multi_system = dquadratization_multi_equilibrium(
        multistable_eq_system, dissipative_equilibrium_list, method='numpy', display=False)
    assert dissipative_multi_system[0] == 67108864 or dissipative_multi_system[0] == 134217728


def test_algorithm_3_large_sympy():
    dissipative_multi_system = dquadratization_multi_equilibrium(
        multistable_eq_system, dissipative_equilibrium_list, method='sympy-naive', display=False)
    assert dissipative_multi_system[0] == 67108864 or dissipative_multi_system[0] == 134217728


def test_algorithm_3_large_Routh():
    dissipative_multi_system = dquadratization_multi_equilibrium(
        multistable_eq_system, dissipative_equilibrium_list, method='Routh-Hurwitz', display=False)
    assert dissipative_multi_system[0] == 67108864 or dissipative_multi_system[0] == 134217728
