import sympy as sp

from functools import reduce
from typing import List, Iterable
import AST_walk as ast
from SymbolsHolder import SymbolsHolder, make_derivative_symbol
from sympy import init_printing
from sympy import symbols
import Combinations as comb

import sympy as sp


class EquationSystem:
    def __init__(self, system):
        self._righthand = []
        self._variables = set()
        self._constants = set()

        self.is_poly = self.check_is_polynomial_system(system)

        if self.is_poly is False:
            raise ValueError("The system is not a polynomial system, please enter a polynomial system.")

        for eq in system:
            self._righthand.append(eq.rhs)
            self._variables.update(eq.lhs.free_symbols)
            self._constants.update(eq.rhs.free_symbols)

        self._constants = self._constants - self._variables
        self.all_terms_RHS = self.get_all_terms_RHS(system)
        self.V = self._variables | {1}

        # build the dictionary of variables and equations
        self._dict_variables_equations = self.get_dict_variables_equations(
            system)
        self.VSquare = self.get_VSquare(system) | self.variables
        self.NSquare = self.all_terms_RHS - (self.VSquare & self.all_terms_RHS)
        self.NQuadratic = comb.get_all_nonquadratic(self.V)
        self.degree = comb.max_degree_monomial(self.all_terms_RHS)
        self.dimension = len(self.variables)

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
    
    def check_is_polynomial_system(self, system):
        for eq in system:
            if ast.is_polynomial_function_all(eq.rhs) is not True:
                return False
        return True

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
        decomposition = term.expand().as_ordered_factors()
        for factor in decomposition:
            if factor in self._constants:
                constant_terms.add(factor)
        return constant_terms
    
    def find_number_coefficient(self, term):
        num_coefficient = 1
        for key, value in term.as_coefficients_dict().items():
            num_coefficient = value
            break
        return num_coefficient

    def get_all_terms_RHS(self, system):
        all_terms_RHS = set()
        for eq in system:
            rhs = eq.rhs
            for term in sp.expand(rhs).as_ordered_terms():
                # Here we use combinations to find all the terms in RHS, we first
                term_dict = term.as_powers_dict()
                list_of_terms = []
                for key in term_dict:
                    if key in self.variables:
                        list_of_terms.append(key**term_dict[key])
                all_terms_RHS.add(reduce(lambda x, y: x * y, list_of_terms))

        return all_terms_RHS

    
