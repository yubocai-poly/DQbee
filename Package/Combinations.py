import sympy as sp
import cmath


def score_function(monomial):
    """
    Role: Compute the score function of a monomial
    Input: monomial without coefficient
    Output: score function of the monomial
    """
    powers = monomial.as_powers_dict()
    num = 1
    for var, exp in powers.items():
        num *= (exp + 1)
    return num


def degree_function(monomial):
    """
    Role: Compute the degree of a monomial
    """
    if monomial.is_constant():
        return 0
    return sum(sp.degree_list(monomial))


def set_to_score_dict(monomial_set):
    """
    Role: Compute the score function of a monomial set
    Input: a set of monomials
    Output: a dictionary with key is monomial and value is score function of the monomial
    """
    score_dict = {}
    for monomial in monomial_set:
        score = score_function(monomial)
        score_dict[monomial] = score
    return score_dict


def get_decompositions(monomial):
    """
    Role: Compute all the decompositions of a monomial
    Input: a tuple (2, 1, 2) represent the degree of a monomial x^2y^1z^2
    Output: a set of decompositions of the monomial, in tuple form
    """
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
    """
    Transform a monomial to tuple form
    Input: a monomial
    Output: tuple form of the monomial
    Example: x^2y^1z^2 -> (2, 1, 2)
    """
    monomial_dict = monomial.as_powers_dict()
    monomial_variables = list(monomial_dict.keys())
    monomial = tuple(monomial_dict.values())
    return monomial_variables, monomial


def tuple_to_monomial(monomial_variables, monomial):
    """
    Transform a tuple to monomial form
    """
    exponents = list(monomial)
    return sp.prod([monomial_variables[i]**exponents[i] for i in range(len(monomial_variables))])


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


def check_quadratic(variables, insert_variable: sp.Poly):
    """"
    This function check if the insert_variable is quadratic in the system
    """
    for variable in variables:
        divisor, remainder = sp.div(insert_variable, variable)
        if remainder == 0:
            if divisor in variables:
                return True
    return False


def get_all_nonquadratic(variables):
    """
    This function returns all the nonquadratic variables in the system
    """
    nonquadratic = set()
    for variable in variables:
        if variable != 1 and degree_function(variable) != 1 and check_quadratic(variables, variable) == False:
            nonquadratic.add(variable)
    return nonquadratic
