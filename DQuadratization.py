import sympy as sp
from EquationSystem import EquationSystem
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
    new_variables_dict.update(new_variables_dict_copy)

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

        
    return ODQ_system, Introduced_variables


    

            
    



             

