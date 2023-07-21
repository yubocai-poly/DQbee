import sympy as sp
from EquationSystem import EquationSystem
from sympy import symbols, diff, expand
from sympy.core.relational import Equality
import copy
from Combinations import *

def calculate_polynomial_derivative(polynomial, equations: EquationSystem):
    derivative = 0
    variables = polynomial.free_symbols & equations.variables

    for variable in variables:
        variable_derivative = diff(polynomial, variable)
        equation_variable_diff = equations.dict_variables_equations[variable]
        derivative += variable_derivative * equation_variable_diff

    return sp.expand(derivative)

def calculate_new_subproblem(equations: EquationSystem, subproblem=None):
    equations = copy.deepcopy(equations)
    new_system = equations.system
    for el in subproblem:
        new_rhs = calculate_polynomial_derivative(el, equations)
        new_system.append(Equality(el, new_rhs))
    return EquationSystem(new_system)

