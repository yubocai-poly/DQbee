import sys
import time
sys.path.append("../../")
from Package.DQuadratization import *
from Package.EquationSystem import *
from Package.Combinations import *
import sympy as sp

# Total number of oscillators
N = 8

# Creating variables for positions and velocities
xs = sp.symbols(" ".join([f"x{i}" for i in range(N)]))
dxs = sp.symbols(" ".join([f"dx{i}" for i in range(N)]))

# Defining the coupling and dumping parameters
def A_entry(i, j):
    if i == j:
        return 1
    if i - j in [-1, 1]:
        return sp.Rational(1, 3)
    return 0

A = sp.Matrix([[A_entry(i, j) for j in range(N)] for i in range(N)])
delta = 2

# Defining the equations
eqs = [sp.Eq(xs[i], dxs[i]) for i in range(N)]
for i in range(N):
    linear = (A * sp.Matrix(xs))[i]
    eqs.append(sp.Eq(dxs[i], linear - linear**3 - delta * dxs[i]))
coupled_duffing = EquationSystem(eqs)

# Computing inner quadratization only
start = time.time()
inner_quadratization = optimal_inner_quadratization(coupled_duffing)
end = time.time()
print(f"Inner quadratization computed in {end - start} seconds")

# Computing the dissipative equilibria
equilibria = []
for rhs in range(2**N):
    rhs_vect = []
    for i in range(N):
        if rhs % 2 == 1:
            rhs_vect.append(1)
        else:
            rhs_vect.append(-1)
        rhs = rhs // 2
    equilibria.append(list(A.solve(sp.Matrix(rhs_vect))) + [0 for _ in range(N)])
    
# Computing the dissipative quadratization
for method in ['numpy']:#['numpy', 'Routh-Hurwitz']:
    start = time.time()
    dissipative_multi_system = dquadratization_multi_equilibrium(
        coupled_duffing, 
        equilibria, 
        method=method
    )
    end = time.time()
    print(f"Diffipative quadratization computed in {end - start} seconds ({method} method)")
