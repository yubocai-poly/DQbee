import sympy as sp
import sys 
sys.path.append("..")
from Package.EquationSystem import EquationSystem


def score_function(polynomial):
    """
    Here the polynomial of input need to be polynomial without the coefficient
    """
    powers = polynomial.as_powers_dict()
    num = 1
    for var, exp in powers.items():
        num *= (exp + 1)
    return num


def find_number_coefficient(term):
    num_coefficient = 1
    for key, value in term.as_coefficients_dict().items():
        num_coefficient = value
        break
    return num_coefficient


def degree_function(polynomial):
    if polynomial.is_constant():
        return 0
    return sum(sp.degree_list(polynomial))


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
            a, b = tuple(list(r[0]) + [i]
                         ), tuple(list(r[1]) + [monomial[-1] - i])
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
        decomposition1 = tuple_to_monomial(
            monomial_variables, decomposition[0])
        decomposition2 = tuple_to_monomial(
            monomial_variables, decomposition[1])
        result.append((decomposition1, decomposition2))
    return result


def max_degree_monomial(dict):
    max_degree = 0
    for monomial in dict:
        degree_monomial = sum(sp.degree_list(monomial))
        if degree_monomial > max_degree:
            max_degree = degree_monomial
    return max_degree


def check_Non_quadratic(variables, insert_variable: sp.Poly):
    """"
    This function check if the insert_variable is quadratic in the system
    """
    for variable in variables:
        divisor, remainder = sp.div(insert_variable, variable)
        if remainder == 0:
            if divisor in variables:
                return False
    return True


def get_all_nonquadratic(variables):
    """
    This function returns all the nonquadratic variables in the system
    """
    nonquadratic = set()
    for variable in variables:
        if variable != 1 and degree_function(variable) != 1 and check_Non_quadratic(variables, variable):
            nonquadratic.add(variable)
    return nonquadratic


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
    system_degree = system.degree
    selected_variable = select_decompose_variable(system)
    decomposition = decomposition_monomial(selected_variable[0])
    valid_decomposition = []

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


def compute_largest_eigenvalue(matrix: sp.Matrix):
    """
    This function computes the largest eigenvalue of a matrix
    """
    eigenvalues = matrix.eigenvals()
    eigenvalues = list(eigenvalues.keys())
    if len(eigenvalues) == 0:
        return 0
    else:
        return max(eigenvalues)
