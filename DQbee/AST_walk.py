"""
Code from Qbee: https://github.com/AndreyBychkov/QBee/tree/a87aabd1fc15717118e4f0f95d0fae3a43ac480e
Author: Andrey Bychkov
Edited by: Yubo Cai
"""

import sympy as sp


def is_polynomial_function(expr: sp.Expr):
    # test one term is polynomial or not
    if expr.is_Add or expr.is_Mul or _is_positive_numeric_power(expr):
        return True
    else:
        return False


def _is_positive_numeric_power(expr: sp.Expr):
    if not expr.is_Pow:
        return False

    power = expr.args[1]
    return power.is_Number and power > 0


def find_non_polynomial(expr: sp.Expr, mode="backward"):
    if mode == "backward":
        return find_non_polynomial_backward(expr)
    if mode == "forward":
        return find_non_polynomial_forward(expr)


def find_non_polynomial_forward(expr: sp.Expr):
    if not is_polynomial_function(expr) and not expr.is_Symbol and not expr.is_Number:
        return expr

    results = map(find_non_polynomial_forward, expr.args)
    results = filter(lambda e: e is not None, results)
    results = list(results)

    if results:
        return results[0]
    return None


def find_non_polynomial_backward(expr: sp.Expr):
    if expr.args:
        results = map(find_non_polynomial_backward, expr.args)
        results = filter(lambda e: e is not None, results)
        results = list(results)
        if results:
            return results[0]

    if not is_polynomial_function(expr) and not expr.is_Symbol and not expr.is_Number:
        return expr
    return None


def is_polynomial_function_all(expr: sp.Expr):
    # test all terms are polynomial or not
    if find_non_polynomial(expr) is None:
        return True
    else:
        return (False, find_non_polynomial(expr))
