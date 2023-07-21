import sympy as sp
from EquationSystem import EquationSystem
from sympy import symbols, diff, expand
from sympy.core.relational import Equality
import copy
from Combinations import *
from DifferientialPoly import *
        
def getOneDQuadratization(system: EquationSystem):
    root = copy.deepcopy(system)
    d = system.degree
    while not (len(root.NQuadratic) == 0 and len(root.NSquare) == 0):
        subproblem = decompose_variable(root, d)[1][0]
        root = calculate_new_subproblem(root, subproblem)
    return root, root.variables


def DissipativeQuadratization(system: EquationSystem, d):
    if (len(system.NSquare) == 0 and len(system.NQuadratic) == 0):
        print(system.variables)
        return 
    
    subproblem_set = decompose_variable(system, d)[1]
    # print(subproblem_set)
    # print(system.variables)
    subproblem_list = []
    for subproblem in subproblem_set:
        subproblem_list.append(calculate_new_subproblem(system, subproblem))
        # print(subproblem_list[-1].variables)
    for subproblem in subproblem_list:
        DissipativeQuadratization(subproblem, d)


            
    



             

