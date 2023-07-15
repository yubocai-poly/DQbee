import sympy as sp
from EquationSystem import EquationSystem
from sympy import symbols, diff, expand
from sympy.core.relational import Equality

def calculate_polynomial_derivative(polynomial, equations: EquationSystem):
    derivative = 0
    variables = polynomial.free_symbols & equations.variables

    for variable in variables:
        variable_derivative = diff(polynomial, variable)
        equation_variable_diff = equations.dict_variables_equations[variable]
        derivative += variable_derivative * equation_variable_diff

    return sp.expand(derivative)