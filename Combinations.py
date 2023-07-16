import sympy as sp
from sympy import symbols
from EquationSystem import EquationSystem


def score_function(polynomial):
    """
    Here the polynomial of input need to be polynomial without the coefficient
    """
    powers = polynomial.as_powers_dict()
    num = 1
    for var, exp in powers.items():
        num *= (exp + 1)
    return num


def set_to_score_dict(polynomial_set):
    score_dict = {}
    for polynomial in polynomial_set:
        score = score_function(polynomial)
        score_dict[polynomial] = score
    return score_dict


def get_decompositions(monomial):
    if len(monomial) == 0:
        return {(tuple(), tuple())}
    result = set()
    prev_result = get_decompositions(tuple(monomial[:-1]))
    for r in prev_result:
        for i in range(abs(monomial[-1]) + 1):
            i = i if monomial[-1] >= 0 else -i
            a, b = tuple(list(r[0]) + [i]), tuple(list(r[1]) + [monomial[-1] - i])
            result.add((min(a, b), max(a, b)))
    return result

def monomial_to_tuple(monomial):
    monomial_dict = monomial.as_powers_dict()
    monomial_variables = list(monomial_dict.keys())
    monomial = tuple(monomial_dict.values())
    return monomial_variables, monomial

def tuple_to_monomial(monomial_variables, monomial):
    monomial_dict = {}
    for i in range(len(monomial_variables)):
        monomial_dict[monomial_variables[i]] = monomial[i]
    monomial = sp.Mul(*[key**value for key, value in monomial_dict.items()])
    return monomial

def decomposition_monomial(monomial):
    result = []
    monomial_variables, monomial = monomial_to_tuple(monomial)
    decompositions = get_decompositions(monomial)
    for decomposition in decompositions:
        decomposition1 = tuple_to_monomial(monomial_variables, decomposition[0])
        decomposition2 = tuple_to_monomial(monomial_variables, decomposition[1])
        result.append((decomposition1, decomposition2))
    return result
    

def max_degree_monomial(dict):
    max_degree = 0
    for monomial in dict:
        degree_monomial = sum(sp.degree_list(monomial))
        if degree_monomial > max_degree:
            max_degree = degree_monomial
    return max_degree

    
def check_Non_quadratic(system: EquationSystem, insert_variable: sp.Poly):
    """"
    This function check if the insert_variable is quadratic in the system
    """
    V = system.variables
    for el in V:
        divisor, remindor = sp.div(insert_variable, el)
        if remindor == 0:
            if divisor in V:
                return True
    return False