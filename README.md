# Dissipative-Quadratiation-Package

## 1. Introduction to the Problem

<kbd>Definition 1 (Quadratization).</kbd> For a system (1), by a quadratization is a pair consisting of
- a list of new variables

$$
y_1=g_1(\mathbf{x}), \ldots, y_m=g_m(\mathbf{x})
$$

- and two lists

$$
\mathbf{q}_1(\mathbf{x}, \mathbf{y})=(q_{1,1}(\mathbf{x}), \ldots, q_{1, n}(\mathbf{y})) \quad \text { and } \quad \mathbf{q}_2(\mathbf{x}, \mathbf{y})=(q_{2,1}(\mathbf{x}, \mathbf{y}), \ldots, q_{2, m}(\mathbf{x}, \mathbf{y}))
$$

of $m+n$-variate polynomials in $\mathbf{x}$ and $\mathbf{y}=\left(y_1, \ldots, y_m\right)$ such that the degree of each of of $\mathbf{q}_1$ and $\mathbf{q}_2$ is at most two and

$$
\mathbf{x}^{\prime}=\mathbf{q}_1(\mathbf{x}, \mathbf{y}) \quad \text { and } \quad \mathbf{y}^{\prime}=\mathbf{q}_2(\mathbf{x}, \mathbf{y}) .
$$

If all the polynomials $g_1, \ldots, g_m$ are monomials, the quadratization is called monomial quadratization.

<kbd>Definition 2 (Equilibrium).</kbd> For a polynomial ODE system (1), a point $\mathbf{x}^* \in$ $\mathbb{R}^n$ is called an equilibrium if $\mathbf{p}\left(\mathbf{x}^*\right)=0$.

<kbd>Definition 3 (Dissipativity).</kbd> An ODE system (1) is called dissipative at an equilibrium point $\mathbf{x}^*$ if all the eigenvalues of the Jacobian $J(\mathbf{p})|_{\mathbf{x}=\mathbf{x}^*}$ of $\mathbf{p}$ and $\mathbf{x}^*$ have negative real part.

<kbd>Definition 4 (Dissipative quadratization).</kbd> Assume that a system (1) is dissipative at an equilibrium point $\mathbf{x}^* \in \mathbb{R}^n$. Then a quadratization given by $\mathbf{g}, \mathbf{q}_1$ and $\mathbf{q}_2$ (see Definition [) is called dissipative at $\mathbf{x}^*$ if the system

$$
\mathbf{x}^{\prime}=\mathbf{q}_1(\mathbf{x}, \mathbf{y}), \quad \mathbf{y}^{\prime}=\mathbf{q}_2(\mathbf{x}, \mathbf{y})
$$

is dissipative at a point $\left(\mathbf{x}^*, \mathbf{g}\left(\mathbf{x}^*\right)\right)$.

## 2. Package Description