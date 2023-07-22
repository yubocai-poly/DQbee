import sympy as sp
import sys 
sys.path.append("..")
from functools import reduce
from typing import List, Iterable
import Package.AST_walk as ast
import Package.Combinations as comb
from IPython.display import Latex

import sympy as sp


class EquationSystem:
    def __init__(self, system):
        self._system = system
        self._righthand = []
        self._variables = set()
        self._constants = set()

        self.is_poly = self.check_is_polynomial_system(system)

        if self.is_poly is False:
            raise ValueError("The system is not a polynomial system, please enter a polynomial system.")

        for eq in system:
            self._righthand.append(eq.rhs)
            self._variables.add(eq.lhs)
            # add the variables like x1 ** 3 as well
            # if comb.degree_function(eq.lhs) > 1:
            #     new_variables_dict = eq.lhs.as_powers_dict()
            #     # combine the new_variables_dict together and add the new variables
            #     list_of_terms = []
            #     for key in new_variables_dict:
            #         list_of_terms.append(key ** new_variables_dict[key])
            #     print(list_of_terms)
            #     self._variables.add(reduce(lambda x, y: x * y, list_of_terms))
            # else:
            #     self._variables.update(eq.lhs.free_symbols)
            self._constants.update(eq.rhs.free_symbols)

        self._constants = self._constants - self._variables
        self.all_terms_RHS = self.get_all_terms_RHS(system)
        self.V = self._variables | {1}

        # build the dictionary of variables and equations
        self._dict_variables_equations = self.get_dict_variables_equations(
            system)
        self.VSquare = self.get_VSquare(system) | self.variables
        self.NSquare = self.all_terms_RHS - (self.VSquare & self.all_terms_RHS)
        self.NQuadratic = comb.get_all_nonquadratic(self.variables)
        self.degree = comb.max_degree_monomial(self.all_terms_RHS)
        self.dimension = len(self.variables)


    @property
    def system(self):
        return self._system
    
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

    def print_system_latex(self):
        latex_system = []
        system = self.system
        for eq in system:
            lhs = sp.latex(eq.lhs)
            rhs = sp.latex(eq.rhs)
            latex_eq = f"({lhs})' &= {rhs}\\\\" if eq != system[-1] else f"({lhs})' &= {rhs}"
            latex_system.append(latex_eq)

        return '\n'.join(latex_system)

    def show_system_latex(self):
        latex_str = f"\\begin{{cases}}\n{self.print_system_latex()}\n\\end{{cases}}"
        return Latex(rf"$${latex_str}$$")

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

    # def find_constant_coefficient(self, term):
    #     constant_terms = set()
    #     decomposition = term.expand().as_ordered_factors()
    #     for factor in decomposition:
    #         if factor in self._constants:
    #             constant_terms.add(factor)
    #     return constant_terms


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
    


    
