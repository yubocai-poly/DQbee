import sympy as sp

from functools import reduce
from typing import List, Iterable
from AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder, make_derivative_symbol
from sympy import init_printing


# class EquationSystem:
#     def __init__(self, equations: List[sp.Eq],
#                  parameter_variables: Iterable[sp.Symbol] = None,
#                  input_variables: Iterable[sp.Symbol] = None):
#         self._equations = equations.copy()
#         self._original_equation_indexes = list(range(len(equations)))
#         self._replacement_equations = list()

#         self._righthand = [eq.args[1] for eq in equations]
#         self._variables = [eq.args[0] for eq in equations]
#         self._constant = list()

#         for eq in self._righthand:
#             # TODO

#     @property
#     def equations(self):
#         return self._equations

#     @property
#     def righthand(self):
#         return self._righthand

#     @property
#     def variables(self):
#         return self._variables

#     @property
#     def constant(self):
#         return self._constant

import sympy as sp

class EquationSystem:
    def __init__(self, system):
        self._righthand = []
        self._variables = set()
        self._constants = set()

        for eq in system:
            self._righthand.append(eq.rhs)
            self._variables.update(eq.lhs.free_symbols)
            self._constants.update(eq.rhs.free_symbols)

        self._constants = self._constants - self._variables

    @property
    def righthand(self):
        return self._righthand

    @property
    def variables(self):
        return self._variables

    @property
    def constants(self):
        return self._constants

    def print_system_latex(self, system):
        latex_system = []

        for eq in system:
            lhs = sp.latex(eq.lhs)
            rhs = sp.latex(eq.rhs)
            latex_eq = f"{lhs}'={rhs}"
            latex_system.append(latex_eq)

        return '\n'.join(latex_system)