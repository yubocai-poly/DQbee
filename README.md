# Dissipative-Quadratiation-Package

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
We will demonstrate usage of Package on the example below. Other interactive examples you can find more in our [example notebook]().

### 3.1. Importing the Package
