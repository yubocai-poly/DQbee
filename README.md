# DQbee

`DQbee` is a python package for computing the dissipative quadratization of polynomial ODE system. 

### Installation
```bash
pip install DQbee
```

## 1. Introduction to the Problem

<kbd>Definition 1 (Quadratization).</kbd> Consider a system of ODEs $\mathbf{x}'=\mathbf{p}(\mathbf{x})$

$$
\begin{cases} 
x_1^{\prime}=f_1(\bar{x})\\
\ldots \\
x_n^{\prime}=f_n(\bar{x})
\end{cases}
\tag{1}
$$

where $\bar{x}=\left(x_1, \ldots, x_n\right)$ and $f_1, \ldots, f_n \in \mathbb{C}[\mathbf{x}]$. Then a list of new variables

$$
y_1=g_1(\bar{x}), \ldots, y_m=g_m(\bar{x})
$$

is said to be a quadratization of (1) if there exist polynomials $h_1, \ldots, h_{m+n} \in$ $\mathbb{C}[\bar{x}, \bar{y}]$ of degree at most two such that
- $x_i^{\prime}=h_i(\bar{x}, \bar{y})$ for every $1 \leqslant i \leqslant n$ which we define as $\mathbf{x}^{\prime}=\mathbf{q}_1(\mathbf{x}, \mathbf{y})$
- $y_j^{\prime}=h_{j+n}(\bar{x}, \bar{y})$ for every $1 \leqslant j \leqslant m$ which we define as $\mathbf{y}^{\prime}=\mathbf{q}_2(\mathbf{x}, \mathbf{y})$

Here we call the number $m$ as the **order of quadratization**. The **optimal quadratization** is the approach that produce the smallest possible order.

<kbd>Definition 2 (Equilibrium).</kbd> For a polynomial ODE system (1), a point $\mathbf{x}^* \in$ $\mathbb{R}^n$ is called an equilibrium if $\mathbf{p}\left(\mathbf{x}^*\right)=0$.

<kbd>Definition 3 (Dissipativity).</kbd> An ODE system (1) is called dissipative at an equilibrium point $\mathbf{x}^\ast$ if all the eigenvalues of the Jacobian $J(\mathbf{p})|_{\mathbf{x}=\mathbf{x}^\ast}$ of $\mathbf{p}$ and $\mathbf{x}^\ast$ have negative real part.

<kbd>Definition 4 (Dissipative quadratization).</kbd> Assume that a system (1) is dissipative at an equilibrium point $\mathbf{x}^* \in \mathbb{R}^n$. Then a quadratization given by $\mathbf{g}, \mathbf{q}_1$ and $\mathbf{q}_2$ is called dissipative at $\mathbf{x}^*$ if the system

$$
\mathbf{x}^{\prime}=\mathbf{q}_1(\mathbf{x}, \mathbf{y}), \quad \mathbf{y}^{\prime}=\mathbf{q}_2(\mathbf{x}, \mathbf{y})
$$

is dissipative at a point $\left(\mathbf{x}^\ast, \mathbf{g}\left(\mathbf{x}^\ast\right)\right)$.

## 2. Package Description
This python package is mainly used to compute the inner-quadratic quadratization and dissipative quadratization of multivariate high-dimensional polynomial ODE system. This algorithm also has tools to perform reachability analysis of the system. Here are some of the specific features of this package:

- Compute the inner-quadratic quadratization of polynomial ODE system.
- Compute the dissipative quadratization of polynomial ODE system at origin.
- Compute the dissipative quadratization of polynomial ODE system at an equilibrium point (input by the user), or at multiple equilibrium points (input by the user) with method `numpy`, `sympy`, or `Routh-Huritwz criterion`.
- Compute the weakly nonlinearity bound $\left| \mathcal{X_{0}} \right|$ of dissipative system.

## 3. Package Usage
We will demonstrate usage of Package on the example below. Other interactive examples you can find more in our [[example file]](https://github.com/yubocai-poly/DQbee/tree/main/Examples). **We recommend to run the package in IPython environment or Jupyter notebook in order to see the result display with $\LaTeX$**.

### 3.1. Importing the Package
```python
import sympy as sp
from DQbee.EquationSystem import *
from DQbee.DQuadratization import *
```

### 3.2. Polynomial ODE System Setting
We take the following duffing equation as an example:

$$
x'' = kx + ax^{3} + bx'
$$

which described damped oscillator with non-linear restoring force. The equation can be written as a first-order system by introducing $x_1 := x, x_2 := x'$ as follows

$$
\begin{cases}
    x_1' = x_2,\\
    x_2' = kx_1 + a x_1^3 + b x_2.
\end{cases}
$$

We can define the system as follows:
```python
x1, x2 = sp.symbols('x1 x2')
k, a, b = sp.symbols('k a b')

system = [
    sp.Eq(x1, x2),
    sp.Eq(x2, k*x1 + a*x1**3 + b*x2)
]

eq_system = EquationSystem(system) # This step transform the system into `EquationSystem` object
display_or_print(eq_system.show_system_latex()) # Display the system in latex format, inner method of `EquationSystem` object
```

### 3.3. Quadratization
In our package, we provide several quadratization functions to the `EquationSystem` object:
- `optimal_inner_quadratization`: compute the inner-quadratic quadratization of a given polynomial ODE system. *(Implementation of Algorithm 1 in the paper)*
- `optimal_dissipative_quadratization`: transform an inner-quadratic system into dissipative quadratization at origin, where the user can choose the value of the coefficients of the linear terms in order to make the system dissipative.
- `dquadratization_one_equilibrium`: compute the dissipative quadratization of a given polynomial ODE system at a given equilibrium point.
- `dquadratization_multi_equilibrium`: compute the dissipative quadratization of a given polynomial ODE system at multiple equilibrium points. *(Implementation of Algorithm 2 in the paper)*

```py
inner_result = optimal_inner_quadratization(eq_system,
                                            display=False) # Set `display` to True if you want to display the result
# Output:
oiq_system = inner_result[0]
sub_oiq_system = inner_result[1]
monomial_to_quadratic_form = inner_result[2] # This is only used for optimal_dissipative_quadratization` function
map_variables = inner_result[3]
```
Here, we will explain with these outputs. We show transformation from `oiq_system` to `sub_oiq_system` by substitution of variables:

$$
\begin{cases}
\left(x_1\right)^{\prime}=x_2 \\
\left(x_2\right)^{\prime}=a x_1^3+b x_2+k x_1 \\
\left(x_1^2\right)^{\prime}=2 x_1 x_2
\end{cases}
\xRightarrow{x_1^2 = w_1} 
\begin{cases}
\left(x_1\right)^{\prime} & =x_2 \\ 
\left(x_2\right)^{\prime} & =a w_1 x_1+b x_2+k x_1 \\ 
\left(w_1\right)^{\prime} & =2 x_1 x_2
\end{cases}
$$

where `map_variables` is the map from the original variables to the new variables. In this case, we have $w_{1}=x_{1}^{2}$. With all the information above, we can compute the dissipative quadratization of the system at origin as follows:
```py
dissipative_result = optimal_dissipative_quadratization(eq_system,
                                                        sub_oiq_system,
                                                        monomial_to_quadratic_form,
                                                        map_variables,
                                                        Input_diagonal=None, # the coefficients of the linear terms
                                                        display=False) 
```
The function will transfer all the linear terms in the ODE of introduced variables into quadratic terms, and keep negative linear terms in order to have a dissipative system *(by keep all diagonal terms of the matrix associated to linear part negative )*. If the `input_diagonal` is None, then the function will choose the largest eigenvalue of orginal system as the diagonal value. In our case, the output is as follows:

$$
\begin{cases}
\left(x_1\right)^{\prime} & =x_2 \\ 
\left(x_2\right)^{\prime} & =a w_1 x_1+b x_2+k x_1 \\ 
\left(w_1\right)^{\prime} & =w_1\left(\frac{b}{2}-\frac{\sqrt{b^2+4 k}}{2}\right)-x_1^2\left(\frac{b}{2}-\frac{\sqrt{b^2+4 k}}{2}\right)+2 x_1 x_2
\end{cases}
\xRightarrow{\text{matrix associated to linear part}} 
\left[\begin{array}{ccc}
0 & 1 & 0 \\
k & b & 0 \\
0 & 0 & \frac{b}{2}-\frac{\sqrt{b^2+4 k}}{2}
\end{array}\right]
$$

We take another example to show how to compute the dissipative quadratization with multiple given equilibrium point. We take the following system as an example:

$$
x' = -x (x - 1) (x - 2)
$$

where $a$ is positive. The system has three equilibrium points: $x_1 = 0, x_2 = 1, x_3 = 2$ but only dissipative at $x_1 = 0, 2$. We can define the system as follows:
```py
multistable_eq_system = EquationSystem([sp.Eq(x, - x * (x - 1) * (x - 2))])
```
Then we can compute the dissipative quadratization of the system at multiple equilibrium points as follows:
```py
result = dquadratization_multi_equilibrium(multistable_eq_system, 
                                           equilibrium = [[0],[2]],
                                           display=True,
                                           method='numpy' # or 'sympy' or 'Routh-Hurwitz'
                                        )
```
which will return the $\lambda$ and the dissipative quadratization of the system which make all the equilibrium points dissipative. In this case, the output is as follows:

$$
\left\{\begin{array}{l}
(x)^{\prime}=-g_1 x+3 x^2-2 x \\
\left(g_1\right)^{\prime}=-2 g_1^2+6 g_1 x-8 g_1+4 x^2
\end{array}\right.
$$

In this function, we provide three methods to compute the dissipative quadratization of the system at multiple equilibrium points:
- `numpy`: use `numpy.linalg.eigvals` to compute the eigenvalues of the Jacobian matrix of the system at equilibrium points, which is the fastest method.
- `sympy`: use `sympy.Matrix.eigenvals` to compute the eigenvalues of the Jacobian matrix of the system at equilibrium points.
- `Routh-Hurwitz`: use Routh-Hurwitz criterion to justify the stability of the system at equilibrium points.

## 4. Artifact Evaluation
We provide the artifact evaluation for the paper. We list the code resourses and the corresponding results in the following list:
- Figure 1 in Example 2 in the paper: [[Jupyter Notebook]](https://github.com/yubocai-poly/DQbee/blob/main/Examples/simulation_unstable.ipynb).

- Implementation of Algorithm 1 in the paper: [[code]](https://github.com/yubocai-poly/DQbee/blob/main/DQbee/DQuadratization.py#L88).

- Implementation of Algorithm 2 in the paper: [[code]](https://github.com/yubocai-poly/DQbee/blob/main/DQbee/DQuadratization.py#L591).

- Example 5 in the paper: [[Jupyter Notebook]](https://github.com/yubocai-poly/DQbee/blob/main/Examples/MultiStability/multistability_system_2.ipynb).

- Section 5.1 - Application to reachability analysis: [[Julia Code]](https://github.com/yubocai-poly/DQbee/tree/main/Examples/Reachability). Plots of the results are shown in the file [[results]](https://github.com/yubocai-poly/DQbee/tree/main/Examples/Reachability/results). In order to execute the program of the reachability analysis, please run the [[start.ipynb]](https://github.com/yubocai-poly/DQbee/blob/main/Examples/Reachability/startup.ipynb).

- Section 5.2 - Small bistable model: resource is shown in the file [[Bistability]](https://github.com/yubocai-poly/DQbee/tree/main/Examples/BiStability).
- Section 5.3 - Coupled Duffing oscillators: [[Python Code]](https://github.com/yubocai-poly/DQbee/blob/main/Examples/CoupledDuffing/coupled_duffing.py).