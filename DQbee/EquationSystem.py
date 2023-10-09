import sys
sys.path.append('../')
from collections import Counter
from IPython.display import Latex
import DQbee.Combinations as comb
from DQbee.Combinations import *
import DQbee.AST_walk as ast
from typing import List, Iterable
from functools import reduce
import copy
import sympy as sp
from sympy import symbols, diff, expand
from sympy.core.relational import Equality


class EquationSystem:
    def __init__(self, system):
        self._system = system
        self._variables = set()
        self._constants = set()

        self.is_poly = self.check_is_polynomial_system(system)

        if self.is_poly is False:
            raise ValueError(
                "The system is not a polynomial system, please enter a polynomial system.")

        for eq in system:
            self._variables.add(eq.lhs)
            self._constants.update(eq.rhs.free_symbols)

        # here constants are the variables that are not in the LHS of the system
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
        self.pruning_rule_OQ_num = self.pruning_rule_OQ()
        self.pruning_rule_ODQ_num = self.pruning_rule_ODQ()

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

    # ---------------------- Methods ----------------------

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
            rhs = sp.latex(eq.rhs.as_expr())
            latex_eq = f"({lhs})' &= {rhs}\\\\" if eq != system[-1] else f"({lhs})' &= {rhs}"
            latex_system.append(latex_eq)

        return '\n'.join(latex_system)

    def max_degree_monomial(self, monomials_set):
        """
        This function compute the degree of the highest degree monomial in the system
        """
        max_degree = 0
        for monomial in monomials_set:
            degree_monomial = comb.degree_function(monomial)
            if degree_monomial > max_degree:
                max_degree = degree_monomial
        return max_degree

    def show_system_latex(self):
        latex_str = f"\\begin{{cases}}\n{self.print_system_latex()}\n\\end{{cases}}"
        return Latex(rf"$${latex_str}$$")

    def get_dict_variables_equations(self):
        return {eq.lhs: eq.rhs for eq in self.system}

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
        """
        all_terms_RHS = set()
        variables = list(self.variables)
        for eq in system:
            rhs = eq.rhs
            terms_dict = rhs.as_poly(variables).as_dict()
            for term_tuple, coef in terms_dict.items():
                # combine the term together
                # Yubo: Here all the terms are in monomial form
                term = sp.prod(
                    [var**exp for var, exp in zip(variables, term_tuple)])
                all_terms_RHS.add(term)

        return all_terms_RHS

    # ---------------------- Pruning Rules ----------------------

    def computation_D_mult(self):
        NS = self.NSquare
        V = self.V
        multiset = []
        for m in NS:
            for v in V:
                divisor, reminder = sp.div(m, v)
                if reminder == 0:
                    multiset.append(divisor)
        return multiset, Counter(multiset).most_common()

    def pruning_rule_OQ(self):
        """
        algorithms 2 in optimal quadratization paper
        """
        k = 0
        len_NS = len(self.NSquare)
        D_mult = self.computation_D_mult()[1]
        while len_NS > sum([i[1] for i in D_mult[:k]]) + k * (k + 1) / 2:
            k += 1
        return k

    def computation_H_mult(self):
        NQ = self.NQuadratic
        V = self.V
        multiset = []
        for m in NQ:
            for v in V:
                if m != v:
                    divisor, reminder = sp.div(m, v)
                    if reminder == 0:
                        multiset.append(divisor)
        return multiset, Counter(multiset).most_common()

    def pruning_rule_ODQ(self):
        """
        algorithm 2 for optimal dissipative quadratization
        """
        k = 0
        len_NQ = len(self.NQuadratic)
        len_NS = len(self.NSquare)
        # print(len_NQ+len_NS)
        D = self.computation_D_mult()[0]
        H = self.computation_H_mult()[0]
        F = H + D
        F_mult = Counter(F).most_common()
        while len_NQ + len_NS > sum([i[1] for i in F_mult[:k]]) + k * (k + 1) / 2:
            k += 1
        return k

    # ---------------------- Quadratization Part ----------------------

    def select_decompose_variable(self):
        """
        This function selects the variable with the lowest score function from NS and NQ to decompose. We also need to check which set of the variable comes from
        """
        NS = self.NSquare
        NQ = self.NQuadratic
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

    def process_decomposition(self, decomposition, d):
        V = self.V
        valid_decomposition = []
        for decompose in decomposition:
            res = []
            for i in range(2):
                if degree_function(decompose[i]) < d and decompose[i] not in V:
                    res.append(decompose[i])
            if len(res) == 0:
                continue
            elif len(res) == 2 and res[0] == res[1]:
                valid_decomposition.append([res[0]])
            elif len(res) == 2 and res[0] != res[1]:
                valid_decomposition.append(res)
            else:
                valid_decomposition.append(res)
        return valid_decomposition

    def decompose_variable(self, d):
        """
        Role: Find all the valid decompositions of the selected variable
        """
        selected_variable = self.select_decompose_variable()
        decomposition = decomposition_monomial(selected_variable[0])
        valid_decomposition = []

        if selected_variable[1] == 'NS':
            # Filter out variables with degree >= system_degree
            valid_decomposition = self.process_decomposition(
                decomposition, d)

        elif selected_variable[1] == 'NQ':
            # Filter out variables with degree >= system_degree
            decomposition = [
                decompose for decompose in decomposition if 1 not in decompose]
            valid_decomposition = self.process_decomposition(
                decomposition, d)

        return valid_decomposition

    def calculate_monomial_derivative(self, polynomial):
        derivative = 0
        variables = polynomial.free_symbols & self.variables
        for variable in variables:
            variable_derivative = diff(polynomial, variable)
            equation_variable_diff = self.dict_variables_equations[variable]
            derivative += variable_derivative * equation_variable_diff
        return derivative.as_poly(list(self.variables))

    def update_system_aux(self, new_variable):
        """
        This function create a new system with the adding of a new variable
        Input: old system, adding variable
        Output: new system
        """
        new_rhs = self.calculate_monomial_derivative(new_variable)
        eq = Equality(new_variable, new_rhs)
        new_system = copy.deepcopy(self)
        new_system._variables.add(new_variable)
        new_system.degree = self.degree
        new_system.V = new_system.variables | {1}

        # update the all_terms_RHS and the degree of the system
        terms_dict = new_rhs.as_poly(list(self.variables)).as_dict()
        for term_tuple, coef in terms_dict.items():
            new_term = reduce(
                lambda x, y: x * y, [var**exp for var, exp in zip(self.variables, term_tuple)])
            new_system.all_terms_RHS.add(new_term)
            if sum(term_tuple) > new_system.degree:
                new_system.degree = sum(term_tuple)

        # update other part
        new_system._dict_variables_equations[new_variable] = new_rhs
        for el in new_system.variables:
            new_system.VSquare.add(el * new_variable)
        new_system.VSquare = new_system.VSquare | new_system.variables
        new_system.NSquare = new_system.all_terms_RHS - \
            (new_system.VSquare & new_system.all_terms_RHS)
        new_system.dimension = len(new_system.variables)
        new_system._system = new_system.system + [eq]
        new_system.pruning_rule_OQ_num = new_system.pruning_rule_OQ()
        new_system.pruning_rule_DOQ_num = new_system.pruning_rule_ODQ()

        # update the NQuadratic
        is_inner_quadratic = comb.check_quadratic(
            new_system.variables, new_variable)
        for el in self.NQuadratic:
            if comb.check_quadratic(new_system.variables, el) == False:
                continue
            else:
                new_system.NQuadratic.remove(el)
        if is_inner_quadratic is False:
            new_system.NQuadratic.add(new_variable)

        return new_system

    def update_system(self, subproblem):
        subproblem = list(subproblem)
        system = copy.deepcopy(self)
        if len(subproblem) == 1:
            return self.update_system_aux(subproblem[0])
        else:
            for i in range(len(subproblem)):
                system = system.update_system_aux(subproblem[i])
            return system
