import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ..EquationSystem import *
from ..DQuadratization import *
from typing import Optional, List, Dict
from sympy import lambdify


def input_information(system: EquationSystem):
    """
    This function gives the information of input variables of a system of equations for simulation.
    """
    print('------------------------------------ Input of coefficients ------------------------------------')
    if system.constants is None:
        print('The system only contains numerical constants.')
    else:
        print('The system contains the following symbolic constants:')
        print('Please give the value of the constants in the form of a dictionary, otherwise the default value will be given.')
        print(system.constants)

    print('-------------------------- Input of initial state of the variables ---------------------------')
    print('Please give the initial state of the variables in the form of a dictionary, otherwise the default value will be given.')
    print(system.variables)
    return


def system_to_odeint(system: EquationSystem,
                     initial_state: Optional[Dict] = None,
                     symbolic_args: Optional[Dict] = None,
                     t: Optional[List] = None,
                     map_variables: Optional[Dict] = None):
    """
    This function is used to convert a system of equations to a function that can be used by odeint.
    """
    # Input information
    input_information(system)

    # building the value for symbolic constants and initial state
    if symbolic_args is None:
        symbolic_args = {}
    if initial_state is None:
        initial_state = {}

    for symbol in system.constants:
        if symbol not in symbolic_args.keys():
            symbolic_args[symbol] = 1

    for symbol in system.variables:
        if symbol not in initial_state.keys():
            initial_state[symbol] = 0.1

    if map_variables is None:
        pass
    else:
        for intro_var, expr in map_variables.items():
            # {w1: x1 ** 2} means w1 = x1 ** 2, so I need to compute the initial state of w1 based on the initial state of x1
            if expr.is_Symbol:
                initial_state[intro_var] = initial_state[expr]
            else:
                initial_state[intro_var] = float(expr.subs(initial_state))

    if t is None:
        t = np.linspace(0, 10, 1000)
    else:
        t = np.linspace(t[0], t[1], t[2])

    # Convert the equations to a function that can be used by odeint
    equations = [eq.rhs for eq in system.system]
    f = lambdify([system.variables, system.constants], equations, "numpy")

    # Convert the initial state and symbolic arguments to arrays in the correct order
    y0_with_index = np.array([(var, initial_state[var])
                             for var in system.variables])
    variable_index = [y0_with_index[i][0]
                      for i in range(len(y0_with_index))]
    y0 = np.array([y0_with_index[i][1] for i in range(len(y0_with_index))])
    args = np.array([symbolic_args[const] for const in system.constants])

    # Define a function that can be used by odeint
    def func(y, t):
        return f(y, args)

    # Solve the system using odeint
    solution = odeint(func, y0, t)

    return solution, variable_index, t
