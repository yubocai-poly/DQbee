import sys
sys.path.append('../')
from tbcontrol.symbolic import routh
import copy
import numpy as np
import sympy as sp
from DQbee.EquationSystem import *
from DQbee.Combinations import *
from functools import reduce
from sympy import symbols, diff, expand
from sympy.core.relational import Equality
from sympy import symbols, total_degree
from IPython.display import Latex


def display_or_print(data):
    """
    find out which environment the code is running in and display or print the data
    """
    if 'ipykernel_launcher.py' in sys.argv[0]:
        from IPython.display import display
        display(data)
    else:
        print(data)


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


def inner_quadratization_aux(system: EquationSystem, d, current_optimal):
    """
    Role: This function finds the optimal innter quadratization of a system
    Input: system, dimension of the original system
    Output: the optimal inner quadratization of the system
    """
    optimal_result = current_optimal

    if len(system.NSquare) == 0 and len(system.NQuadratic) == 0:
        if optimal_result == None:
            optimal_result = system.variables
        elif len(system.variables) < len(optimal_result):
            optimal_result = system.variables
        return optimal_result

    subproblem_set = system.decompose_variable(d)
    subproblem_list = []
    for subproblem in subproblem_set:
        pruning_lower_bound = len(
            system.variables) + max(system.pruning_rule_ODQ_num, system.pruning_rule_OQ_num)
        naive_lower_bound = len(subproblem) + len(system.variables)
        if (optimal_result is not None) and (naive_lower_bound >= len(optimal_result) or pruning_lower_bound >= len(optimal_result)):
            continue
        subproblem_list.append(system.update_system(subproblem))
        # subproblem_list.append(calculate_new_subproblem(system, subproblem))
    for subproblem in subproblem_list:
        optimal_result = inner_quadratization_aux(
            subproblem, d, optimal_result)

    return optimal_result


def inner_quadratization(system: EquationSystem, d):
    return inner_quadratization_aux(system, d, None)


def optimal_inner_quadratization(system: EquationSystem, display=True):
    """
    Role: find the optimal inner quadratization of a system
    Input: 
        - system
    Output: 
        - oiq_system: the optimal inner quadratization of the system
        - sub_oiq_system: the optimal system after substituting the new introduced variables
        - monomial_to_quadratic_form: the sub dictionary from monomial_2_quadra which only contain the quadratic expression of new introduced variables, this output is used for Optimal DQ function
        - map_variables: the dictionary of the mapping between the new introduced variables and the original variables, e.g. {w1: x1 ** 2}
    """
    d = system.degree
    optimal_result = inner_quadratization(system, d)
    # Here monomial_2_quadra is the map from monomial to it's quadratic form, this includes the monomial expression of introduced variables to the new introduced variables, eg. x1 ** 2 -> w1, and includes the monomial expression to it quadratic form, eg. x1 ** 2 * x2 -> w1 * x2
    monomial_2_quadra = {}
    monomial_to_quadratic_form = {}
    map_variables = {}
    substitute_system = []
    oiq_variables = optimal_result
    introduced_variables = oiq_variables - system.variables
    oiq_system = system.update_system(introduced_variables)
    if display:
        print('The Original System is: ')
        display_or_print(system.show_system_latex())
        print('The Optimal Dissipative Quadratization is: ')
        display_or_print(oiq_system.show_system_latex())

    # for each new introduced variable, we create a new symbol corresponding to it, like w_1, w_2, ...
    num = 0
    for variable in oiq_variables:
        if variable in system.variables:
            monomial_2_quadra[variable] = variable
            num += 1
        else:
            monomial_2_quadra[variable] = symbols(
                'w' + str(len(monomial_2_quadra) + 1 - num))
            map_variables[monomial_2_quadra[variable]] = variable

    # print the new introduced variables in latex
    if display:
        print('The new introduced variables are: ')
    new_variables_latex = ''
    for variable in introduced_variables:
        new_variables_latex = new_variables_latex + \
            f"{sp.latex(monomial_2_quadra[variable])} = {sp.latex(variable)} \\\\ "
    latex_ = f"\\begin{{cases}}\n{new_variables_latex}\n\\end{{cases}}"
    if display:
        display_or_print(Latex(rf"$${latex_}$$"))

    # Here monomial_to_quadratic_form which is the map from monomial to it's quadratic form only
    # monomial_2_quadra also contains the variable to itself
    monomial_2_quadra_copy = copy.deepcopy(monomial_2_quadra)
    for variable1 in oiq_variables:
        for variable2 in oiq_variables:
            if variable1 in introduced_variables or variable2 in introduced_variables:
                monomial_2_quadra[variable1 * variable2] = monomial_2_quadra_copy[variable1] * \
                    monomial_2_quadra_copy[variable2]
                monomial_to_quadratic_form[
                    variable1 * variable2] = monomial_2_quadra_copy[variable1] * monomial_2_quadra_copy[variable2]
    # uncomment the following line if not make the whole rhs quadratic
    # monomial_2_quadra.update(monomial_2_quadra)

    variables = list(oiq_system.variables)
    for equation in oiq_system.system:
        lhs = equation.lhs
        rhs = equation.rhs.as_poly(variables)
        new_lhs = None
        new_rhs = 0
        # substitute the left hand side to make the degree is 1
        if lhs in monomial_2_quadra_copy:
            new_lhs = monomial_2_quadra_copy[lhs]
        # now we substitute the right handside to make sure it quadratic
        for key, value in rhs.as_dict().items():
            cleaned_term = reduce(
                lambda x, y: x*y, [var ** exp for var, exp in zip(variables, key)])
            if degree_function(cleaned_term) > 2:
                new_rhs += value * \
                    monomial_2_quadra[cleaned_term]
            else:
                new_rhs += value * cleaned_term

        # Here we give the output of the substitution, that's why I use sp.Expr to make sure the output is a sympy expression
        # Yubo: Also, the substitution for sp.Poly is not stable, it occurs an error of not substituting the once when I run the code for ten times
        substitute_system.append(Equality(new_lhs, new_rhs))

    sub_oiq_system = EquationSystem(substitute_system)
    if display:
        print('The Optimal Quadratic Dissipative System is (with substitution): ')
        display_or_print(sub_oiq_system.show_system_latex())

    return oiq_system, sub_oiq_system, monomial_to_quadratic_form, map_variables

# ------------------ Dissipative Quadratization ------------------


def optimal_dissipative_quadratization(original_system: EquationSystem,
                                       iq_system: EquationSystem,
                                       variables_dict: dict,
                                       map_variables: dict,
                                       Input_diagonal=None,
                                       display=True):
    """
    Role: Compute the Optimal Dissipative Quadratization of a system
    Input:
        - original_system: the original system
        - iq—system: the optimal inner-quadratic system from the original system, need to be the one with substitution
        - variables_dict: the dictionary of the new introduced variables, e.g. {x1 ** 2: w1}
        - map_variables: the dictionary of the mapping between the new introduced variables and the original variables, e.g. {w1: x1 ** 2}
        - Input_diagonal: the list of values of the diagonal of F1 that we want to input, if None, we just use the largest eigenvalue of F1_original as default
    Output:
        - the optimal dissipative quadratization of the system
        - Latex print of the system
        - Latex print of F1
    """
    if not (len(iq_system.NQuadratic) == 0 and len(iq_system.NSquare) == 0):
        raise ValueError(
            'The system (second input) is not a inner-quadratic, please enter a inner-quadratic system.')
    n_original = original_system.dimension
    n = iq_system.dimension
    F1_original = sp.zeros(n_original)
    F1 = sp.zeros(n)
    list_term = []
    type_system = None
    original_variables = list(original_system.variables)
    for equation in iq_system.system:
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

    if display:
        print('------------------------------ Optimal Dissipative Quadratization ------------------------------')
    if original_system.constants != set():
        if display:
            print('The system contains symbolic constants. Therefore, we cannot make sure that the original system is dissipative.')
            print('If the system can not jusity which eigenvalue is the largest one due to the symbolic constants, the program will choose the first one as the default value for the largest eigenvalue.')
        type_system = 'symbolic'
    else:
        type_system = 'numeric'

    largest_eigenvalue = compute_largest_eigenvalue(F1_original, type_system)
    if display:
        print('The eigenvalue with the largest real part is (real part): ',
              largest_eigenvalue)

    if Input_diagonal == None:
        Input_diagonal = [largest_eigenvalue] * n
        if display:
            print(
                'Mode: Default diagonal value, we will choose the largest eigenvalue as the diagonal value.')
    elif len(Input_diagonal) != n - n_original:
        raise ValueError(
            f'The length of the input diagonal is {len(Input_diagonal)}, but the length of the diagonal should be {n - n_original}')
    elif max(Input_diagonal) >= 0:
        raise ValueError(
            'The diagonal value should be all negative in order to make the system dissipative.')

    F1[0:n_original, 0:n_original] = F1_original
    for i in range(n_original, n):
        F1[i, i] = Input_diagonal[i - n_original]

    new_system = iq_system.system[0:n_original]
    for i in range(n_original, n):
        equation = iq_system.system[i]
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
    if display:
        print('The converted Optimal Dissipative Quadratization System is: ')
        display_or_print(dissipative_system.show_system_latex())
        print('The matrix  associated to the linear part system is:')
        display_or_print(Latex(rf"$${sp.latex(F1)}$$"))

    return dissipative_system, F1, map_variables

# -------------------------------------------------------------------
# Compute the weakly nonlinearity bound for |x|
# -------------------------------------------------------------------


def compute_weakly_nonlinearity(system: EquationSystem, display=True):
    """
    Role: Compute the bound |x| for a system being weakly nonlinear
    """
    if display:
        print("-------------------------- Compute Weakly Nonlinearity bound for |x| --------------------------")
        print("-------------------------- Warning: Please enter a quadratic ODE system --------------------------")
    n = system.dimension
    F1 = sp.zeros(n, n)
    F2 = sp.zeros(n, (n ** 2))
    list_term = []
    upper_bound = 0
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
    if display:
        display_or_print(Latex(rf"$$F_{1}={sp.latex(F1)}$$"))
        display_or_print(Latex(rf"$$F_{2}={sp.latex(F2)}$$"))

    largest_eigen_F1 = np.abs(
        max([complex(eigenvalue).real for eigenvalue in F1.eigenvals().keys()]))
    operator_norm_F2 = np.abs(F2.norm(2))
    expression = r"The system is said to be weakly nonlinear if the ratio $$R:=\frac{\left\|X_0\right\|\left\|F_2\right\|}{\left|\Re\left(\lambda_1\right)\right|} < 1 \Rightarrow 0 < \|X_0\| < \frac{\left|\Re\left(\lambda_1\right)\right|}{\left\|F_2\right\|}$$"
    if display:
        display_or_print(Latex(expression))
    if largest_eigen_F1 == 0:
        if display:
            display_or_print(Latex(
                rf"The bound for $\|X\|$ is invalid since $\left|\Re\left(\lambda_1\right)\right|=0$"))
    else:
        upper_bound = operator_norm_F2 / largest_eigen_F1
        if display:
            display_or_print(Latex(
                rf"The upper bound for $\|X\|$ is: $$0 < \|X\| < {upper_bound}$$ where $\left|\Re\left(\lambda_1\right)\right|={largest_eigen_F1}$ and $\left\|F_2\right\|={operator_norm_F2}$"))
    return F2, upper_bound

# -------------------------------------------------------------------
# Algorithm 2 - Compute the dissipative quadratization with one equilibrium
# -------------------------------------------------------------------


def system_transformation(system: list,
                          equilibrium: dict):
    """
    Role: This function perform the coordinate transformation on the system at a given equilibrium point
    Input:
        - system: the input system, represented as a list of equations
        - equilibrium: the equilibrium point, need to indicate the value of all the variables
    Output:
        - transformed_system: the system after coordinate transformation
        - coordinate_transformation: the dictionary of the coordinate transformation, e.g. {'u1': y - 3}
    """
    # Create a substitution dictionary where each variable is replaced with a new symbol
    substitution = {var: sp.Symbol(f'u{index}') for index, var in enumerate(
        equilibrium.keys(), start=1)}

    # Create a coordinate transformation dictionary, e.g., {'u1': y - 3}
    coordinate_transformation = {var: new_var + value for new_var,
                                 (var, value) in zip(substitution.values(), equilibrium.items())}

    # Create a new system where each equation has been coordinate-transformed
    transformed_system = []
    for eq in system:
        # Only replace variables on the left side
        lhs = eq.lhs.subs(substitution)
        # Replace variables and constants on the right side
        rhs = eq.rhs.subs(coordinate_transformation).expand()
        transformed_system.append(sp.Eq(lhs, rhs))

    return transformed_system, {v: k - equilibrium[k] for k, v in substitution.items()}


def dquadratization_one_equilibrium(system: EquationSystem,
                                    equilibrium: dict,
                                    display=True):
    """
        Role: The system compute the dissipative quadratization of a system at a given equilibrium point (Algorithm 2 in the paper)
        Input:
            - system: the input system, represented as a list of equations
            - equilibrium: the equilibrium point, need to indicate the value of all the variables
        Output:
            - odq_transformed_eq_system: the dissipative quadratization of the system after coordinate transformation
            - coordinate_transformation: the dictionary of the coordinate transformation, e.g. {'u1': y - 3}
            - oiq_transformed_eq_system: the optimal inner-quadratic system of the system after coordinate transformation
    """
    # checking the input, wheter the number of equilibrium is equal to the number of variables or not
    if len(equilibrium) != system.dimension:
        raise ValueError(
            'The dimension of the equilibrium point is not equal to the number of variables in the system, please check the input')

    if display:
        print("-------------------------- Compute a quadratization dissipative at a given equilibrium --------------------------")
        print("The system before coordinate transformation is: ")
        display_or_print(system.show_system_latex())
        print("-------------------------- Dissipative Quadratization on transformed system --------------------------")
    transformed_system, coordinate_transformation = system_transformation(
        system.system, equilibrium)
    transformed_eq_system = EquationSystem(transformed_system)
    oiq_transformed_eq_system = optimal_inner_quadratization(
        transformed_eq_system)
    odq_transformed_eq_system = optimal_dissipative_quadratization(
        transformed_eq_system, oiq_transformed_eq_system[1], oiq_transformed_eq_system[2], oiq_transformed_eq_system[3])

    return odq_transformed_eq_system, coordinate_transformation, oiq_transformed_eq_system

# -------------------------------------------------------------------
# Algorithm 3 - Compute the dissipative quadratization with multiple equilibrium
# -------------------------------------------------------------------


def sort_equation_system(system: EquationSystem,
                         map_variables: dict):
    """
    Role: sort the equation system based on the degree of the variables
    Input:
        - system: the input system
        - map_variables: a dictionary of the new introduced variables to their expression, e.g. {w1: x1 ** 2}
    """
    sorted_system = [0] * system.dimension
    degree_dict = {k: total_degree(v.as_poly())
                   for k, v in map_variables.items()}
    sorted_variables = sorted(degree_dict, key=degree_dict.get)
    new_variables = symbols(
        ','.join(f'g{i+1}' for i in range(len(sorted_variables))))
    if not isinstance(new_variables, tuple):
        new_variables = [new_variables]
    new_variables_expression = {
        new_variables[i]: map_variables[sorted_variables[i]] for i in range(len(sorted_variables))}
    origin_to_new_variables = {sorted_variables[i]: new_variables[i] for i in range(
        len(sorted_variables))}  # from old to new

    index = 0
    for equation in system.system:
        lhs = equation.lhs
        rhs = equation.rhs
        if lhs not in map_variables.keys():
            rhs = rhs.subs(origin_to_new_variables).expand()
            sorted_system[index] = sp.Eq(lhs, rhs)
            index += 1
        else:
            new_lhs = origin_to_new_variables[lhs]
            new_rhs = system.dict_variables_equations[lhs].subs(
                origin_to_new_variables).expand()
            position = new_variables.index(new_lhs) + index
            sorted_system[position] = sp.Eq(new_lhs, new_rhs)

    new_eq_system = EquationSystem(sorted_system)
    return new_eq_system, new_variables_expression


def innerquadratic_representation(map_variables: dict):
    """
    Role: Since the input map_variables is inner-quadratic, therefore, we try to find the quadratic form of the variables, e.g. {g1: x1 ** 2, g2: x1 ** 3} => {g1: x1 ** 2, g2: x1 * g1}
    """
    variable_to_quadratic_form = {}
    map_variables_inverse = {v: k for k, v in map_variables.items()}
    # the degree of variables is increasing
    map_variables_value_list = list(map_variables.values())

    for variable, expression in map_variables.items():
        if degree_function(expression) == 2:
            variable_to_quadratic_form[variable] = expression
        else:
            for el in map_variables_value_list:
                divisor, reminder = sp.div(expression, el)
                if reminder == 0 and divisor != 1:
                    if degree_function(divisor) == 1:
                        variable_to_quadratic_form[variable] = divisor * \
                            map_variables_inverse[el]
                        break
                    elif divisor in map_variables_value_list:
                        variable_to_quadratic_form[variable] = map_variables_inverse[divisor] * \
                            map_variables_inverse[el]
                        break
                else:
                    continue

    return variable_to_quadratic_form


def equilibrium_list_to_dict(system: EquationSystem,
                             equilibrium: list,
                             map_variables: dict):
    """
    Role: This function match the list of equilibrium points to the dictionary of equilibrium points in order to use sympy.subs for computation, eg. [[0, 0], [1, 1]] for [x1, x2] => {x1: 0, x2: 0} and {x1: 1, x2: 1}
    """
    equilibrium_dict = []
    n = system.dimension
    list_lhs = [lhs for lhs in system.dict_variables_equations.keys()]
    for equilibria in equilibrium:
        equilibria_dict = {}
        if len(equilibria) != n:
            raise ValueError(
                'The dimension of the equilibrium point is not equal to the number of variables in the system, please check the input')
        for i in range(n):
            equilibria_dict[list_lhs[i]] = equilibria[i]
        # add the new introduced variables to the equilibrium point
        for key, value in map_variables.items():
            equilibria_dict[key] = value.subs(equilibria_dict)
        equilibrium_dict.append(equilibria_dict)
    return equilibrium_dict


def aux_sympy_naive(list_jacobian_subs_equilirbium, lambda_):
    lambda_value = 0
    while True:
        # plug in the value
        all_eigenvalues_negative = True
        for jacobian_subs_equilibrium in list_jacobian_subs_equilirbium:
            jacobian_matrix_value = jacobian_subs_equilibrium.subs(
                lambda_, lambda_value)
            eigenvalues = list(jacobian_matrix_value.eigenvals().keys())
            max_real_eigen = max(
                [complex(eigenvalue).real for eigenvalue in eigenvalues])
            if max_real_eigen >= 0:
                # if the largest real part eigenvalue is not negative, then we need to increase the lambda value
                all_eigenvalues_negative = False
                if lambda_value == 0:
                    lambda_value = 1
                else:
                    lambda_value *= 2
                break
        if all_eigenvalues_negative:
            break

    return lambda_value


def aux_numpy(list_jacobian_subs_equilirbium, lambda_):
    # finish the computation with numpy
    # first transform the sympy matrix to numpy matrix, contain the lambda symbol
    # plug in the lambda value
    # Create a numpy function for each jacobian matrix
    numpy_jacobians = [sp.lambdify(lambda_, jacobian_matrix_value, 'numpy')
                       for jacobian_matrix_value in list_jacobian_subs_equilirbium]

    lambda_value = 0
    while True:
        # Plug in the lambda value
        all_eigenvalues_negative = True
        for numpy_jacobian in numpy_jacobians:
            jacobian_matrix_value = numpy_jacobian(lambda_value)
            eigenvalues = np.linalg.eigvals(jacobian_matrix_value)
            max_real_eigen = max(eigenvalues.real)
            if max_real_eigen >= 0:
                # If the largest real part eigenvalue is not negative, then we need to increase the lambda value
                all_eigenvalues_negative = False
                if lambda_value == 0:
                    lambda_value = 1
                else:
                    lambda_value *= 2
                break
        if all_eigenvalues_negative:
            break

    return lambda_value


def aux_routh_hurwitz(list_jacobian_subs_equilirbium, lambda_):
    # finish the computation with Routh-Hurwitz criterion
    # first transform the sympy matrix into characteristic polynomial
    character_poly_list = [sp.Poly(jacobian_matrix_value.charpoly(
        lambda_), lambda_) for jacobian_matrix_value in list_jacobian_subs_equilirbium]

    # Convert sympy Polynomials to numpy
    routh_arrays = [routh(p) for p in character_poly_list]

    lambda_value = 0
    while True:
        # Plug in the lambda value
        all_eigenvalues_negative = True
        for routh_matrix in routh_arrays:
            # plug in the lambda value and check whether the left hand column must have entries with all the same signs
            _lambda = routh_matrix.free_symbols.pop()
            routh_matrix_lambda = routh_matrix.subs(_lambda, lambda_value)
            if not (np.all(np.sign(routh_matrix_lambda[:, 0]) > 0) or np.all(np.sign(routh_matrix_lambda[:, 0]) < 0)):
                # If the largest real part eigenvalue is not negative, then we need to increase the lambda value
                all_eigenvalues_negative = False
                if lambda_value == 0:
                    lambda_value = 1
                else:
                    lambda_value *= 2
                break
        if all_eigenvalues_negative:
            break

    return lambda_value


def dquadratization_multi_equilibrium(system: EquationSystem,
                                      equilibrium: list,
                                      display=True,
                                      method='numpy'):
    """
    Algorithm 3 in the paper
    Role: The algorithm take a input system and a list of equilibrium points, and compute the dissipative quadratization of the system with the equilibrium points
    Input:
        - system: the input system
        - equilibrium: a list of equilibrium points, each equilibrium point is a list of values of the variables, e.g. [[0, 0], [1, 1]] for [x1, x2] => {x1: 0, x2: 0} and {x1: 1, x2: 1}
        method: The default method is using numpy for computation, or with routh-hurwitz criterion and naiv_sympy
    Output:
        - lambda_value: the lambda value of the system
        - substract_eq_system_with_lambda: the system after substracting the lambda value
    The function will display the dissipaitve quadratizaiton system with the lambda value
    """
    lambda_ = sp.symbols('lambda')
    if display:
        print("-------------------------- Compute a quadratization dissipative with multi given equilibrium --------------------------")
        print("The original system is: ")
        display_or_print(system.show_system_latex())

    innerquadratic_system, map_variables = optimal_inner_quadratization(
        system, display=False)[1], optimal_inner_quadratization(system, display=False)[3]
    sorted_system, map_variables_expression = sort_equation_system(
        innerquadratic_system, map_variables)
    list_of_equilibrium_dict = equilibrium_list_to_dict(
        system, equilibrium, map_variables_expression)
    variable_to_quadratic_form = innerquadratic_representation(
        map_variables_expression)

    # compute the system substract with lambda * (lhs - quadratic form)
    substract_system = []
    rhs_for_jacobian = []
    lhs_list = []
    for eq in sorted_system.system:
        lhs = eq.lhs
        rhs = eq.rhs
        lhs_list.append(lhs)
        if lhs not in map_variables_expression:
            substract_system.append(sp.Eq(lhs, rhs))
            rhs_for_jacobian.append(rhs)
        else:
            new_rhs = rhs - lambda_ * (lhs - variable_to_quadratic_form[lhs])
            substract_system.append(sp.Eq(lhs, new_rhs.expand()))
            rhs_for_jacobian.append(new_rhs.expand())

    substract_eq_system = EquationSystem(substract_system)

    if display:
        print("-------------------------- The Dissipative quadratized system --------------------------")
        display_or_print(substract_eq_system.show_system_latex())

    # compute the jacobian matrix of the system
    jacobian_matrix = sp.Matrix(rhs_for_jacobian).jacobian(lhs_list)
    list_jacobian_subs_equilirbium = []
    for equilibrium in list_of_equilibrium_dict:
        list_jacobian_subs_equilirbium.append(
            jacobian_matrix.subs(equilibrium))

    # three computation methods
    if method == 'numpy':
        lambda_value = aux_numpy(list_jacobian_subs_equilirbium, lambda_)
    elif method == 'Routh-Hurwitz':
        lambda_value = aux_routh_hurwitz(
            list_jacobian_subs_equilirbium, lambda_)
    elif method == 'sympy-naive':
        lambda_value = aux_sympy_naive(list_jacobian_subs_equilirbium, lambda_)
    else:
        raise ValueError('Please enter a valid method name')

    substract_system_with_lambda = [sp.Eq(lhs, rhs.subs(
        lambda_, lambda_value)) for lhs, rhs in zip(lhs_list, rhs_for_jacobian)]
    substract_eq_system_with_lambda = EquationSystem(
        substract_system_with_lambda)
    if display:
        print("-------------------------- The lambda value --------------------------")
        print('The lambda value is: ', lambda_value)
        print("-------------------------- The Dissipative quadratized system with lambda value --------------------------")
        display_or_print(substract_eq_system_with_lambda.show_system_latex())
        print("-------------------------- The Jacobian matrix --------------------------")
        display_or_print(Latex(rf"$$J={sp.latex(jacobian_matrix)}$$"))

    return lambda_value, substract_eq_system_with_lambda
