import sympy as sp
import numpy as np
from EquationSystem import EquationSystem
from functools import reduce
from sympy import symbols, diff, expand
from sympy.core.relational import Equality
import copy
from Combinations import *
from DifferientialPoly import *
from IPython.display import Latex
        
def getOneDQuadratization(system: EquationSystem):
    root = copy.deepcopy(system)
    d = system.degree
    while not (len(root.NQuadratic) == 0 and len(root.NSquare) == 0):
        subproblem = decompose_variable(root, d)[1][0]
        root = calculate_new_subproblem(root, subproblem)
    return root, root.variables


def DissipativeQuadratization(system: EquationSystem, d, optimal_result):
    if len(system.NSquare) == 0 and len(system.NQuadratic) == 0:
        if optimal_result[0] is None:
            optimal_result[0] = system.variables
        elif len(system.variables) < len(optimal_result[0]):
            optimal_result[0] = system.variables
        return
    
    subproblem_set = decompose_variable(system, d)[1]
    subproblem_list = []
    for subproblem in subproblem_set:
        if (optimal_result[0] is not None) and len(subproblem) + len(system.variables) >= len(optimal_result[0]):
            continue
        subproblem_list.append(calculate_new_subproblem(system, subproblem))
    for subproblem in subproblem_list:
        DissipativeQuadratization(subproblem, d, optimal_result)


def OptimalDissipativeQuadratization(system: EquationSystem):
    d = system.degree
    optimal_result = [None]  
    DissipativeQuadratization(system, d, optimal_result)
    new_variables_dict = {}
    new_variables_dict_ = {}
    map_variables = {}
    original_variables_dict = {}
    substitute_system = []
    new_variables_latex = ''
    ODQ_variables = optimal_result[0]
    Introduced_variables = ODQ_variables - system.variables
    ODQ_system = calculate_new_subproblem(system, Introduced_variables)
    print('The Original System is: ')
    display(system.show_system_latex())
    print('The Optimal Dissipative Quadratization is: ')
    display(ODQ_system.show_system_latex())
    
    # for each new introduced variable, we create a new symbol corresponding to it, like w_1, w_2, ...
    for variable in Introduced_variables:
        new_variables_dict[variable] = symbols('w' + str(len(new_variables_dict) + 1))
        map_variables[symbols('w' + str(len(new_variables_dict)))] = variable
    
    # print the new introduced variables in latex
    print('The new introduced variables are: ')
    for variable in Introduced_variables:
        new_variables_latex = new_variables_latex + f"{sp.latex(new_variables_dict[variable])} = {sp.latex(variable)} \\\\ "
    latex_ = f"\\begin{{cases}}\n{new_variables_latex}\n\\end{{cases}}"
    display(Latex(rf"$${latex_}$$"))

    for variable in system.variables:
        original_variables_dict[variable] = variable
    new_variables_dict.update(original_variables_dict)
    new_variables_dict_copy = copy.deepcopy(new_variables_dict)

    for variable1 in ODQ_variables:
        for variable2 in ODQ_variables:
            if variable1 in Introduced_variables or variable2 in Introduced_variables:
                new_variables_dict[variable1 * variable2] = new_variables_dict_copy[variable1] * new_variables_dict_copy[variable2]
                new_variables_dict_[variable1 * variable2] = new_variables_dict_copy[variable1] * new_variables_dict_copy[variable2]
    # uncomment the following line if not make the whole rhs quadratic
    # new_variables_dict.update(new_variables_dict_copy)

    for equation in ODQ_system.system:
        lhs = equation.lhs
        rhs = equation.rhs
        new_lhs = None
        new_rhs = None
        # substitute the left hand side to make the degree is 1
        if lhs in new_variables_dict_copy:
            new_lhs = new_variables_dict_copy[lhs]
        # now we substitute the right handside to make sure it quadratic
        for term in sp.expand(rhs).as_ordered_terms():
            term_dict = term.as_powers_dict()
            list_of_terms = []
            for key in term_dict:
                if key in system.variables:
                    list_of_terms.append(key**term_dict[key])
            cleaned_term = sp.Mul(*list_of_terms)
            if degree_function(cleaned_term) > 2:
                rhs = rhs.subs(cleaned_term, new_variables_dict[cleaned_term])
        new_rhs = rhs
        substitute_system.append(Equality(new_lhs, new_rhs))

    sub_ODQ_system = EquationSystem(substitute_system)
    print('The Optimal Quadratic Dissipative System is (with substitution): ')
    display(sub_ODQ_system.show_system_latex())

        
    return ODQ_system, new_variables_dict, sub_ODQ_system, new_variables_dict_, map_variables

# ------------------ Dissipative Quadratization ------------------

def dissipative_numeric(original_system: EquationSystem,
                        DQ_system: EquationSystem,
                        variables_dict: dict,
                        map_variables: dict):
    n_original = original_system.dimension
    n = DQ_system.dimension
    F1_original = sp.zeros(n_original)
    F1 = sp.zeros(n)
    list_term = []
    for equation in DQ_system.system:
        list_term.append(equation.lhs)
    index = 0
    for equation in original_system.system:
        rhs = equation.rhs
        for terms in sp.expand(rhs).as_ordered_terms():
            if degree_function(terms) == 1:
                # add the coefficient to the matrix, with the position of the term
                coefficient = find_number_coefficient(terms)
                F1_original[index, list_term.index(terms / coefficient)] = coefficient
        index += 1

    largest_eigen = compute_largest_eigenvalue(F1_original)
    if largest_eigen >= 0:
        raise ValueError('The original system is not dissipative.')
    
    F1[0:n_original, 0:n_original] = F1_original
    # for the diagonal part, we need to add the largest eigenvalue to the diagonal
    for i in range(n_original, n):
        F1[i, i] = largest_eigen

    # now we can build the new system
    new_system = DQ_system.system[0:n_original]
    for i in range(n_original, n):
        equation = DQ_system.system[i]
        lhs = equation.lhs
        if map_variables[lhs] not in variables_dict:
            rhs = largest_eigen * lhs + equation.rhs - largest_eigen * map_variables[lhs]
        else:
            rhs = largest_eigen * lhs + equation.rhs - largest_eigen * variables_dict[map_variables[lhs]]
        new_system.append(Equality(lhs, rhs))

    return F1_original, F1, largest_eigen, EquationSystem(new_system)


def dissipative_symbolic(original_system: EquationSystem,
                        DQ_system: EquationSystem,
                        variables_dict: dict,
                        map_variables: dict):
    n_original = original_system.dimension
    n = DQ_system.dimension
    F1_original = sp.zeros(n_original)
    F1 = sp.zeros(n)
    list_term = []
    constants = original_system.constants
    for equation in DQ_system.system:
        list_term.append(equation.lhs)
    index = 0
    for equation in original_system.system:
        rhs = equation.rhs
        for terms in sp.expand(rhs).as_ordered_terms():
            term_dict = terms.as_powers_dict()
            list_of_coefficients = []
            for key in term_dict:
                if key not in original_system.variables:
                    list_of_coefficients.append(key**term_dict[key])
            coefficient = reduce(lambda x, y: x * y, list_of_coefficients) if list_of_coefficients != [] else 1
            term = terms / coefficient
            if degree_function(term) == 1:
                F1_original[index, list_term.index(term)] = coefficient
        index += 1

    largest_eigen = compute_largest_eigenvalue(F1_original)
    if largest_eigen >= 0:
        print('The largest eigenvalue is: ', largest_eigen)
        print('Is the original system dissipative? ', largest_eigen < 0)
        raise ValueError('The original system is not dissipative.')
    else:
        print('The largest eigenvalue is: ', largest_eigen)

    F1[0:n_original, 0:n_original] = F1_original
    # for the diagonal part, we need to add the largest eigenvalue to the diagonal
    for i in range(n_original, n):
        F1[i, i] = largest_eigen

    new_system = DQ_system.system[0:n_original]
    for i in range(n_original, n):
        equation = DQ_system.system[i]
        lhs = equation.lhs
        if map_variables[lhs] not in variables_dict:
            rhs = largest_eigen * lhs + equation.rhs - largest_eigen * map_variables[lhs]
        else:
            rhs = largest_eigen * lhs + equation.rhs - largest_eigen * variables_dict[map_variables[lhs]]
        new_system.append(Equality(lhs, rhs))

    return F1_original, F1, largest_eigen, EquationSystem(new_system)

def make_system_dissipative(original_system: EquationSystem,
                            DQ_system: EquationSystem,
                            variables_dict: dict,
                            map_variables: dict):
    if not (len(DQ_system.NQuadratic) == 0 and len(DQ_system.NSquare) == 0):
        raise ValueError('The system (second input) is not a inner-quadratic, please enter a inner-quadratic system.')
    if original_system.constants != set():
        print('The system contains symbolic constants. Therefore, we cannot make sure that the original system is dissipative \n\n You may see THIS TYPE OF ERROR: "cannot determine truth value of Relational", this means the program can not tell which one is the largest eigenvalue since there are symbolic constants.')
        result = dissipative_symbolic(original_system, DQ_system, variables_dict, map_variables)
        print('The converted Optimal Dissipative Quadratization System is: ')
        display(result[3].show_system_latex())
        print('The matrix  associated to the linear part system is:')
        display(Latex(rf"$$F_{1}={sp.latex(result[1])}$$"))
    else:
        # This system is purely numeric
        result = dissipative_numeric(original_system, DQ_system, variables_dict, map_variables)
        print('The converted Optimal Dissipative Quadratization System is: ')
        display(result[3].show_system_latex())
        print('The matrix  associated to the linear part system is:')
        display(Latex(rf"$$F_{1}={sp.latex(result[1])}$$"))

    return result[3]
        




            
    



             

