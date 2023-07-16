import sympy as sp


def polynomial_to_dict(polynomial):
    if polynomial.is_number:
        return

    if polynomial.is_Symbol:
        return {polynomial: 1}

    if len(polynomial.free_symbols) == 1:
        var = polynomial
        if isinstance(var, sp.Pow):
            base, exp = var.as_base_exp()
            if isinstance(base, sp.Symbol):
                return {base: exp}

    terms = polynomial.expand().as_ordered_terms()
    powers = {}

    for term in terms:
        if isinstance(term, sp.Mul):
            variables = []
            for factor in term.args:
                if isinstance(factor, sp.Pow):
                    base, exp = factor.as_base_exp()
                    if isinstance(base, sp.Symbol):
                        variables.append((base, exp))
                elif isinstance(factor, sp.Symbol):
                    variables.append((factor, 1))

            for var, exp in variables:
                if var in powers:
                    powers[var] += exp
                else:
                    powers[var] = exp

    return powers


def score_function(polynomial):
    """
    Here the polynomial of input need to be polynomial without the coefficient
    """
    powers = polynomial_to_dict(polynomial)
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
