import copy
import numpy as np
import sympy as sp
from Package.EquationSystem import *
from Package.Combinations import *
from functools import reduce
from sympy import symbols, diff, expand
from sympy.core.relational import Equality
from IPython.display import Latex
import sys
sys.path.append("..")
# from Package.DifferientialPoly import *


def compute_largest_eigenvalue(matrix: sp.Matrix, type_system):
    """
    This function computes the largest eigenvalue of a matrix
    """
    eigenvalues = matrix.eigenvals()
    eigenvalues = list(eigenvalues.keys())
    if len(eigenvalues) == 0:
        return 0
    if type_system == 'symbolic':
        try:
            return max([eigenvalue.real for eigenvalue in eigenvalues])
        except:
            return eigenvalues[0]
    else:
        # get the eigenvalue with the largest real part
        max_real_eigen = max(
            [complex(eigenvalue).real for eigenvalue in eigenvalues])
        if max_real_eigen >= 0:
            print(eigenvalues)
            raise ValueError(
                'The largest eigenvalue is not negative, the original system is not dissipative')
        else:
            return max_real_eigen


def inner_quadratization(system: EquationSystem, d, optimal_result):
    """
    Role: This function finds the optimal innter quadratization of a system
    Input: system, dimension of the original system, optimal_result which is a list
    Output: None, but the optimal_result will be updated
    """
    if len(system.NSquare) == 0 and len(system.NQuadratic) == 0:
        if optimal_result[0] is None:
            optimal_result[0] = system.variables
        elif len(system.variables) < len(optimal_result[0]):
            optimal_result[0] = system.variables
        return

    subproblem_set = system.decompose_variable(d)
    subproblem_list = []
    for subproblem in subproblem_set:
        if (optimal_result[0] is not None) and (len(subproblem) + len(system.variables) >= len(optimal_result[0]) or len(system.variables) + max(system.pruning_rule_ODQ_num, system.pruning_rule_OQ_num) >= len(optimal_result[0])):
            continue
        subproblem_list.append(system.update_system(subproblem))
    for subproblem in subproblem_list:
        inner_quadratization(subproblem, d, optimal_result)

# def inner_quadratization(system: EquationSystem, d, current_optimal):
#     """
#     Role: This function finds the optimal innter quadratization of a system
#     Input: system, dimension of the original system
#     Output: the optimal inner quadratization of the system
#     """
#     optimal_result = current_optimal

#     if len(system.NSquare) == 0 and len(system.NQuadratic) == 0:
#         if sys.dimension < current_optimal
#             return system
#          return current_optimal

#     subproblem_set = system.decompose_variable(d)
#     subproblem_list = []
#     for subproblem in subproblem_set:
#         if (optimal_result[0] is not None) and (len(subproblem) + len(system.variables) >= len(optimal_result[0]) or len(system.variables) + max(system.pruning_rule_ODQ_num, system.pruning_rule_OQ_num) >= len(optimal_result[0])):
#             continue
#         subproblem_list.append(system.update_system(subproblem))
#         # subproblem_list.append(calculate_new_subproblem(system, subproblem))
#     for subproblem in subproblem_list:
#         current_optimal = inner_quadratization(subproblem, d, current_optimal)
#         if optimal_result[0] is None or len(new_optimal_result[0]) < len(optimal_result[0]):
#             optimal_result = new_optimal_result

#     return optimal_result


def optimal_inner_quadratization(system: EquationSystem):
    """
    Role: find the optimal inner quadratization of a system
    Input: 
        - system
    Output: 
        - OIQ_system: the optimal inner quadratization of the system
        - sub_OIQ_system: the optimal system after substituting the new introduced variables
        - monomial_to_quadratic_form: the sub dictionary from monomial_2_quadra which only contain the quadratic expression of new introduced variables, this output is used for Optimal DQ function
        - map_variables: the dictionary of the mapping between the new introduced variables and the original variables, e.g. {w1: x1 ** 2}
    """
    d = system.degree
    optimal_result = [None]
    inner_quadratization(system, d, optimal_result)
    # optimal_result = inner_quadratization(system, d)
    monomial_2_quadra = {}
    monomial_to_quadratic_form = {}
    map_variables = {}
    substitute_system = []
    OIQ_variables = optimal_result[0]
    introduced_variables = OIQ_variables - system.variables
    OIQ_system = system.update_system(introduced_variables)
    print('The Original System is: ')
    display(system.show_system_latex())
    print('The Optimal Dissipative Quadratization is: ')
    display(OIQ_system.show_system_latex())

    # for each new introduced variable, we create a new symbol corresponding to it, like w_1, w_2, ...
    num = 0
    for variable in OIQ_variables:
        if variable in system.variables:
            monomial_2_quadra[variable] = variable
            num += 1
        else:
            monomial_2_quadra[variable] = symbols(
                'w' + str(len(monomial_2_quadra) + 1 - num))
            map_variables[monomial_2_quadra[variable]] = variable


    # print the new introduced variables in latex
    print('The new introduced variables are: ')
    new_variables_latex = ''
    for variable in introduced_variables:
        new_variables_latex = new_variables_latex + \
            f"{sp.latex(monomial_2_quadra[variable])} = {sp.latex(variable)} \\\\ "
    latex_ = f"\\begin{{cases}}\n{new_variables_latex}\n\\end{{cases}}"
    display(Latex(rf"$${latex_}$$"))

    # Here monomial_to_quadratic_form which is the map from monomial to it's quadratic form only
    # monomial_2_quadra also contains the variable to itself
    monomial_2_quadra_copy = copy.deepcopy(monomial_2_quadra)
    for variable1 in OIQ_variables:
        for variable2 in OIQ_variables:
            if variable1 in introduced_variables or variable2 in introduced_variables:
                monomial_2_quadra[variable1 * variable2] = monomial_2_quadra_copy[variable1] * \
                    monomial_2_quadra_copy[variable2]
                monomial_to_quadratic_form[
                    variable1 * variable2] = monomial_2_quadra_copy[variable1] * monomial_2_quadra_copy[variable2]
    # uncomment the following line if not make the whole rhs quadratic
    # monomial_2_quadra.update(monomial_2_quadra)

    variables = list(OIQ_system.variables)
    for equation in OIQ_system.system:
        lhs = equation.lhs
        rhs = equation.rhs.expand()
        new_lhs = None
        new_rhs = None
        # substitute the left hand side to make the degree is 1
        if lhs in monomial_2_quadra_copy:
            new_lhs = monomial_2_quadra_copy[lhs]
        # now we substitute the right handside to make sure it quadratic
        for term in sp.expand(rhs).as_ordered_terms():
            term_dict = term.as_poly(variables).as_dict()
            for key, value in term_dict.items():
                cleaned_term = reduce(
                    lambda x, y: x*y, [var**exp for var, exp in zip(variables, key)])
                if degree_function(cleaned_term) > 2:
                    rhs = rhs.subs(
                        cleaned_term, monomial_2_quadra[cleaned_term])
        new_rhs = rhs
        substitute_system.append(Equality(new_lhs, new_rhs))

    sub_OIQ_system = EquationSystem(substitute_system)
    print('The Optimal Quadratic Dissipative System is (with substitution): ')
    display(sub_OIQ_system.show_system_latex())

    return OIQ_system, sub_OIQ_system, monomial_to_quadratic_form, map_variables

# ------------------ Dissipative Quadratization ------------------


def optimal_dissipative_quadratization(original_system: EquationSystem,
                                       IQ_system: EquationSystem,
                                       variables_dict: dict,
                                       map_variables: dict,
                                       Input_diagonal=None):
    """
    Role: Compute the Optimal Dissipative Quadratization of a system
    Input:
        - original_system: the original system
        - IQ—system: the optimal inner-quadratic system from the original system, need to be the one with substitution
        - variables_dict: the dictionary of the new introduced variables, e.g. {x1 ** 2: w1}
        - map_variables: the dictionary of the mapping between the new introduced variables and the original variables, e.g. {w1: x1 ** 2}
        - Input_diagonal: the list of values of the diagonal of F1 that we want to input, if None, we just use the largest eigenvalue of F1_original as default
    Output:
        - the optimal dissipative quadratization of the system
        - Latex print of the system
        - Latex print of F1
    """
    if not (len(IQ_system.NQuadratic) == 0 and len(IQ_system.NSquare) == 0):
        raise ValueError(
            'The system (second input) is not a inner-quadratic, please enter a inner-quadratic system.')
    n_original = original_system.dimension
    n = IQ_system.dimension
    F1_original = sp.zeros(n_original)
    F1 = sp.zeros(n)
    list_term = []
    type_system = None
    original_variables = list(original_system.variables)
    for equation in IQ_system.system:
        list_term.append(equation.lhs)
    index = 0

    for equation in original_system.system:
        rhs = equation.rhs
        term_dict = rhs.as_poly(original_variables).as_dict()
        for term, coef in term_dict.items():
            if sum(term) == 1:
                new_term = reduce(
                    lambda x, y: x * y, [var**exp for var, exp in zip(original_variables, term)])
                F1_original[index, list_term.index(new_term)] = coef
        index += 1

    print('------------------------------ Optimal Dissipative Quadratization ------------------------------')
    if original_system.constants != set():
        print('The system contains symbolic constants. Therefore, we cannot make sure that the original system is dissipative.')
        print('If the system can not jusity which eigenvalue is the largest one due to the symbolic constants, the program will choose the first one as the default value for the largest eigenvalue.')
        type_system = 'symbolic'
    else:
        type_system = 'numeric'

    largest_eigenvalue = compute_largest_eigenvalue(F1_original, type_system)
    print('The eigenvalue with the largest real part is (real part): ',
          largest_eigenvalue)

    if Input_diagonal == None:
        Input_diagonal = [largest_eigenvalue] * n
        print('Mode: Default diagonal value, we will choose the largest eigenvalue as the diagonal value.')
    elif len(Input_diagonal) != n - n_original:
        raise ValueError(
            f'The length of the input diagonal is {len(Input_diagonal)}, but the length of the diagonal should be {n - n_original}')
    elif max(Input_diagonal) >= 0:
        raise ValueError(
            'The diagonal value should be all negative in order to make the system dissipative.')

    F1[0:n_original, 0:n_original] = F1_original
    for i in range(n_original, n):
        F1[i, i] = Input_diagonal[i - n_original]

    new_system = IQ_system.system[0:n_original]
    for i in range(n_original, n):
        equation = IQ_system.system[i]
        lhs = equation.lhs
        if map_variables[lhs] not in variables_dict:
            rhs = Input_diagonal[i - n_original] * lhs + equation.rhs - \
                Input_diagonal[i - n_original] * map_variables[lhs]
        else:
            rhs = Input_diagonal[i - n_original] * lhs + equation.rhs - \
                Input_diagonal[i - n_original] * \
                variables_dict[map_variables[lhs]]
        new_system.append(Equality(lhs, rhs))

    dissipative_system = EquationSystem(new_system)
    print('The converted Optimal Dissipative Quadratization System is: ')
    display(dissipative_system.show_system_latex())
    print('The matrix  associated to the linear part system is:')
    display(Latex(rf"$${sp.latex(F1)}$$"))

    return dissipative_system, F1


def ComputeWeaklyNonlinearity(system: EquationSystem):
    """
    Role: Compute the bound |x| for a system being weakly nonlinear
    """
    print("-------------------------- Compute Weakly Nonlinearity bound for |x| --------------------------")
    print("-------------------------- Warning: Please enter a quadratic ODE system --------------------------")
    n = system.dimension
    F1 = sp.zeros(n, n)
    F2 = sp.zeros(n, (n ** 2))
    list_term = []
    for equation in system.system:
        list_term.append(equation.lhs)

    index = 0
    for equation in system.system:
        rhs = equation.rhs
        term_dict = rhs.as_poly(list_term).as_dict()
        for term, coef in term_dict.items():
            if sum(term) == 1:
                nonzero_indices = [index for index,
                                   value in enumerate(term) if value != 0]
                F1[index, nonzero_indices[0]] = coef
            if sum(term) == 2:
                nonzero_indices = [index for index,
                                   value in enumerate(term) if value != 0]
                if len(nonzero_indices) == 1:
                    F2[index, nonzero_indices[0] * n + nonzero_indices[0]] = coef
                else:
                    F2[index, nonzero_indices[0] * n + nonzero_indices[1]] = coef

        index += 1
    display(Latex(rf"$$F_{1}={sp.latex(F1)}$$"))
    display(Latex(rf"$$F_{2}={sp.latex(F2)}$$"))

    largest_eigen_F1 = np.abs(
        max([complex(eigenvalue).real for eigenvalue in F1.eigenvals().keys()]))
    operator_norm_F2 = np.abs(F2.norm(2))
    expression = r"The system is said to be weakly nonlinear if the ratio $$R:=\frac{\left\|X_0\right\|\left\|F_2\right\|}{\left|\Re\left(\lambda_1\right)\right|} < 1 \Rightarrow 0 < \|X_0\| < \frac{\left|\Re\left(\lambda_1\right)\right|}{\left\|F_2\right\|}$$"
    display(Latex(expression))
    if largest_eigen_F1 == 0:
        display(Latex(
            rf"The bound for $\|X\|$ is invalid since $\left|\Re\left(\lambda_1\right)\right|=0$"))
    else:
        upper_bound = operator_norm_F2 / largest_eigen_F1
        display(Latex(
            rf"The upper bound for $\|X\|$ is: $$0 < \|X\| < {upper_bound}$$ where $\left|\Re\left(\lambda_1\right)\right|={largest_eigen_F1}$ and $\left\|F_2\right\|={operator_norm_F2}$"))
    return F2
