import sympy as sp
from sympy import symbols, diff, expand
from sympy.core.relational import Equality
import sys 
sys.path.append("..")
import copy
from functools import reduce
from typing import List, Iterable
import Package.AST_walk as ast
from Package.Combinations import *
import Package.Combinations as comb
from IPython.display import Latex

import sympy as sp


class EquationSystem:
    def __init__(self, system):
        self._system = system
        self._variables = set()
        self._constants = set()

        self.is_poly = self.check_is_polynomial_system(system)

        if self.is_poly is False:
            raise ValueError("The system is not a polynomial system, please enter a polynomial system.")

        for eq in system:
            self._variables.add(eq.lhs)
            self._constants.update(eq.rhs.free_symbols)

        self._constants = self._constants - self._variables
        self.all_terms_RHS = self.get_all_terms_RHS(system)
        self.V = self._variables | {1}

        # build the dictionary of variables and equations
        self._dict_variables_equations = self.get_dict_variables_equations()
        self.VSquare = self.get_VSquare(system) | self.variables
        self.NSquare = self.all_terms_RHS - (self.VSquare & self.all_terms_RHS)
        self.NQuadratic = comb.get_all_nonquadratic(self.variables)
        self.degree = self.max_degree_monomial(self.all_terms_RHS)
        self.dimension = len(self.variables)


    @property
    def system(self):
        return self._system

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

    def max_degree_monomial(self, var_set):
        """
        This function compute the degree of the highest degree monomial in the system
        """
        max_degree = 0
        for monomial in var_set:
            degree_monomial = sum(sp.degree_list(monomial))
            if degree_monomial > max_degree:
                max_degree = degree_monomial
        return max_degree

    def show_system_latex(self):
        latex_str = f"\\begin{{cases}}\n{self.print_system_latex()}\n\\end{{cases}}"
        return Latex(rf"$${latex_str}$$")

    def get_dict_variables_equations(self):
        return {eq.lhs : eq.rhs for eq in self.system}

    def get_VSquare(self, system):
        VSquare = set()
        variables = list(self.variables)
        for i in range(len(variables)):
            for j in range(i, len(variables)):
                VSquare.add(variables[i] * variables[j])

        return VSquare

    def get_all_terms_RHS(self, system):
        """
        get all the terms in the RHS of the system
        Yubo: corrected with the suggestion
        """
        all_terms_RHS = set()
        variables = list(self.variables)
        for eq in system:
            rhs = eq.rhs
            terms_dict = rhs.as_poly(variables).as_dict()
            for term, coef in terms_dict.items():
                # combine the term together
                new_term = reduce(lambda x, y: x * y, [var**exp for var, exp in zip(variables, term)])
        
                all_terms_RHS.add(new_term)
        
        return all_terms_RHS
    
def select_decompose_variable(system: EquationSystem):
    """
    This function selects the variable with the lowest score function from NS and NQ to decompose. We also need to check which set of the variable comes from
    """
    NS = system.NSquare
    NQ = system.NQuadratic
    NS_score_dict = set_to_score_dict(NS)
    NQ_score_dict = set_to_score_dict(NQ)
    # choose the one with the lowest score function
    if not NS_score_dict and not NQ_score_dict:
        return None, None
    elif not NS_score_dict:
        return min(NQ_score_dict, key=NQ_score_dict.get), 'NQ'
    elif not NQ_score_dict:
        return min(NS_score_dict, key=NS_score_dict.get), 'NS'
    else:
        if min(NS_score_dict.values()) < min(NQ_score_dict.values()):
            return min(NS_score_dict, key=NS_score_dict.get), 'NS'
        else:
            return min(NQ_score_dict, key=NQ_score_dict.get), 'NQ'


def decompose_variable(system: EquationSystem, d):
    """
    Role: Find all the valid decompositions of the selected variable
    """
    selected_variable = select_decompose_variable(system)
    decomposition = decomposition_monomial(selected_variable[0])
    valid_decomposition = []
    """
    yubo: Here I didn't simply the code since if I merge them together, for every decomposition, I need to check if the decomposition is valid or not, which will increase the time complexity
    """
    if selected_variable[1] == 'NS':
        # Filter out variables with degree >= system_degree
        for decompose in decomposition:
            res = []
            for i in range(2):
                if degree_function(decompose[i]) < d and decompose[i] not in system.V:
                    res.append(decompose[i])
            if len(res) == 0:
                continue
            elif len(res) == 2 and res[0] == res[1]:
                valid_decomposition.append([res[0]])
            else:
                valid_decomposition.append(res)

    if selected_variable[1] == 'NQ':
        # Filter out variables with degree >= system_degree
        for decompose in decomposition:
            res = []
            if 1 in decompose:
                continue
            else:
                for i in range(2):
                    if degree_function(decompose[i]) < d and decompose[i] not in system.V:
                        res.append(decompose[i])
                if len(res) == 0:
                    continue
                elif len(res) == 2 and res[0] == res[1]:
                    valid_decomposition.append([res[0]])
                else:
                    valid_decomposition.append(res)

    return selected_variable, valid_decomposition
    
# Gleb: this could be naturally a method of EquationSystem
def calculate_polynomial_derivative(polynomial, equations: EquationSystem):
    derivative = 0
    variables = polynomial.free_symbols & equations.variables

    for variable in variables:
        variable_derivative = diff(polynomial, variable)
        equation_variable_diff = equations.dict_variables_equations[variable]
        derivative += variable_derivative * equation_variable_diff

    return sp.expand(derivative)
    

def update_system_aux(system: EquationSystem, new_variable):
    """
    This function create a new system with the adding of a new variable
    Input: old system, adding variable
    Output: new system
    """
    new_rhs = calculate_polynomial_derivative(new_variable, system)
    eq = Equality(new_variable, new_rhs)
    new_system = copy.deepcopy(system)
    new_system._variables.add(new_variable)
    new_system.degree = system.degree
    new_system.V = new_system.variables | {1}

    # update the all_terms_RHS
    terms_dict = new_rhs.as_poly(list(system.variables)).as_dict()
    for term, coef in terms_dict.items():
        new_term = reduce(lambda x, y: x * y, [var**exp for var, exp in zip(system.variables, term)])
        new_system.all_terms_RHS.add(new_term)
        if sum(term) > new_system.degree:
            new_system.degree = sum(term)
    
    # update other part
    new_system._dict_variables_equations[new_variable] = new_rhs
    for el in new_system.variables:
        new_system.VSquare.add(el * new_variable)
    new_system.VSquare = new_system.VSquare | new_system.variables 
    new_system.NSquare = new_system.all_terms_RHS - (new_system.VSquare & new_system.all_terms_RHS)
    new_system.dimension = len(new_system.variables)
    new_system._system = new_system.system + [eq]


    # update the NQuadratic
    is_inner_quadratic = comb.check_non_quadratic(new_system.variables, new_variable)
    for el in system.NQuadratic:
        if comb.check_non_quadratic(new_system.variables, el):
            continue
        else:
            new_system.NQuadratic.remove(el)
    if is_inner_quadratic:
        new_system.NQuadratic.add(new_variable)

    return new_system

def update_system(system: EquationSystem, subproblem):
    subproblem = list(subproblem)
    if len(subproblem) == 1:
        return update_system_aux(system, subproblem[0])
    else:
        for i in range(len(subproblem)):
            system = update_system_aux(system, subproblem[i])
        return system