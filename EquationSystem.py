import sympy as sp

from functools import reduce
from typing import List, Iterable
from AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder, make_derivative_symbol
from sympy import init_printing

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
        self.all_terms_RHS = self.get_all_terms_RHS(system)

        # build the dictionary of variables and equations
        self._dict_variables_equations = self.get_dict_variables_equations(
            system)

        self.VSquare = self.get_VSquare(system)

    @property
    def righthand(self):
        return self._righthand

    @property
    def variables(self):
        return self._variables

    @property
    def constants(self):
        return self._constants

    @property
    def dict_variables_equations(self):
        return self._dict_variables_equations

    def print_system_latex(self, system):
        latex_system = []

        for eq in system:
            lhs = sp.latex(eq.lhs)
            rhs = sp.latex(eq.rhs)
            latex_eq = f"{lhs}'={rhs}"
            latex_system.append(latex_eq)

        return '\n'.join(latex_system)

    def get_dict_variables_equations(self, system):
        dict_variables_equations = {}
        for eq in system:
            dict_variables_equations[eq.lhs] = eq.rhs
        return dict_variables_equations

    def get_VSquare(self, system):
        VSquare = set()
        variables = list(self.variables)
        for i in range(len(variables)):
            for j in range(i, len(variables)):
                VSquare.add(variables[i] * variables[j])

        return VSquare
    
    def find_constant_coefficient(self, term):
        constant_terms = set()
        decomposition = term.as_ordered_factors()
        for factor in decomposition:
            if factor in self._constants:
                constant_terms.add(factor)
        return constant_terms

    def get_all_terms_RHS(self, system):
        all_terms_RHS = set()
        for eq in system:
            rhs = eq.rhs
            for term in rhs.args:
                # Here we only want to keep the terms with the constant coefficient, which get off the term in self._constants
                constant_terms = self.find_constant_coefficient(term)
                if constant_terms:
                    all_terms_RHS.add(term / reduce(lambda x, y: x * y, constant_terms))
                else:
                    all_terms_RHS.add(term)

        return all_terms_RHS

