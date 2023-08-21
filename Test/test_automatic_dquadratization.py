import pytest
import math
import sys
import sympy as sp
import time

sys.path.append("..")

from Package.DQuadratization import *
from Package.EquationSystem import *
from Package.Combinations import *
from Package.Simulation.numerical import *


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
OIQ_vanderpol, sub_OIQ_vanderpol, monomial_to_quadratic_form, map_variables = optimal_inner_quadratization(
    eq_system_vanderpol)
DIQ_vanderpol, F1, map_variables = optimal_dissipative_quadratization(
    eq_system_vanderpol, sub_OIQ_vanderpol, monomial_to_quadratic_form, map_variables, [-3])


# ---------------- Testing for Quadratization Part ----------------
def test_input_system():
    assert eq_system_vanderpol.system == system_vanderpol


def test_optimal_inner_quadratization():
    test_OIQ_vanderpol_system = [sp.Eq(y, z), sp.Eq(
        z, mu*z*(1 - y**2) - y), sp.Eq(y**2, 2*y*z)]
    test_sub_OIQ_vanderpol_system = [sp.Eq(y, z), sp.Eq(
        z, -mu*w1*z + mu*z - y), sp.Eq(w1, 2*y*z)]
    assert OIQ_vanderpol.system == test_OIQ_vanderpol_system and sub_OIQ_vanderpol.system == test_sub_OIQ_vanderpol_system


def test_dissipative_inner_quadratization():
    test_F1 = sp.Matrix([[0, 1, 0], [-1, mu, 0], [0, 0, -3]])
    test_DIQ_vanderpol_system = [sp.Eq(y, z), sp.Eq(
        z, -mu*w1*z + mu*z - y), sp.Eq(w1, -3*w1 + 3*y**2 + 2*y*z)]
    assert DIQ_vanderpol.system == test_DIQ_vanderpol_system and F1 == test_F1


# ---------------- Testing for Weakly nonlinearity ----------------
def test_computeWeaklyNonlinearity():
    mu = -1
    system_vanderpol_ = [
        sp.Eq(y, z),
        sp.Eq(z, mu * (1 - y ** 2) * z - y)
    ]
    eq_system_vanderpol_ = EquationSystem(system_vanderpol_)
    test_upper_bound = 2 * math.sqrt(13)
    test_F2 = sp.Matrix([[0, 0, 0, 0, 0, 0, 0, 0, 0], [
        0, 0, 0, 0, 0, 1, 0, 0, 0], [3, 2, 0, 0, 0, 0, 0, 0, 0]])
    OIQ_vanderpol_, sub_OIQ_vanderpol_, monomial_to_quadratic_form_, map_variables_ = optimal_inner_quadratization(
        eq_system_vanderpol_)
    DIQ_vanderpol_, F1_, map_variables_ = optimal_dissipative_quadratization(
        eq_system_vanderpol, sub_OIQ_vanderpol_, monomial_to_quadratic_form_, map_variables_, [-3])
    F2, upper_bound = computeWeaklyNonlinearity(DIQ_vanderpol_)
    assert (test_upper_bound - upper_bound) < 1e-4 and F2 == test_F2


# ---------------- Testing for Simulation ----------------
initial_state = {y: 0.1, z: 0}
symbolic_args = {mu: 1}


def test_simulation_original_system():
    state_original = system_to_odeint(eq_system_vanderpol, t=[
                                      0, 20, 1000], symbolic_args=symbolic_args, initial_state=initial_state)
    assert state_original is not None


def test_simulation_DQ_system():
    state_dissipative = system_to_odeint(DIQ_vanderpol, t=[
                                         0, 20, 1000], symbolic_args=symbolic_args, initial_state=initial_state)
    assert state_dissipative is not None


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
