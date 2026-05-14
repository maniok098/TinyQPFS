# TinyQPFS

Tiny Quadratic Programming From Scratch

The focus is simplicity, determinism, and extremely fast repeated solves for fixed sparse equality-constrained QPs, where only the linear terms are updated online.

## Motivation

For applications where computation time is extremely critical — for example, galvanometer scanner control with sampling times on the order of 10–100 μs — solving a full general-purpose QP can be unnecessarily expensive. In many cases, an equality-constrained QP is already sufficient.

(And if box constraints due to actuator saturation etc. are still absolutely necessary... we can always blame out mechanical engineers — or simply buy better drives.)

There is also an inherent trade-off in MPC itself:

* Faster sampling rates require longer prediction horizons to capture system dynamics adequately.
* At the same time, computational resources per control cycle become increasingly limited.

As the sampling frequency increases, evaluating the MPC law can take longer time than the sampling interval itself, while classical PID controllers still remain feasible and robust.

This project explores a simple question:

> Can a tiny, highly specialized sparse QP solver do its job and possibly outperform PID-style control loops at sampling rates around 100 kHz?


## Basic Idea

The core philosophy of TinyQPFS is:

> Perform all expensive numerical work offline.

Instead of running a full optimization procedure online, the solver directly evaluates the KKT system using precomputed sparse factorizations.

For fixed-structure MPC problems where only the linear terms change online:

* The KKT sparsity pattern remains constant
* Matrix factorization can be reused
* Only the RHS vector needs updating at runtime

In this setting, solving the QP becomes little more than evaluating a sparse linear system.


## Problem Formulation

TinyQPFS focuses on equality-constrained quadratic programs of the form:

$$
\min_x ; \frac{1}{2}x^TQx + q^Tx \quad \text{s.t.} \quad Ax=b
$$

The corresponding KKT system is:
```math
\begin{bmatrix}
Q & A^T \\
A & 0
\end{bmatrix}
\begin{bmatrix}
x \\
\lambda
\end{bmatrix}
=
\begin{bmatrix}
-q \\
b
\end{bmatrix}
```

## Implementation Pipeline

### Offline Initialization

All numerically expensive operations are performed once during initialization:

1. Construct the sparse KKT matrix from the equality-constrained QP formulation
2. Apply matrix scaling

   * e.g. Ruiz equilibration
3. Compute AMD ordering for sparsity-preserving factorization
4. Apply sparse LDLᵀ factorization
5. Search for numerically robust regularization parameters

   * e.g. σ and ρ for quasi-definite regularized KKT systems
   * brute-force tuning during initialization is acceptable
6. Precompute and cache all symbolic factorization data

### Online Solve Phase

At runtime, only lightweight operations are performed:

1. Update linear cost term `q`
2. Update equality RHS `b`
3. Update KKT RHS vector
4. Reuse cached sparse LDLᵀ factorization
5. Solve via forward/backward substitution
6. Apply iterative refinement to improve numerical robustness
7. Warm-start from the previous solution

No matrix rebuilding or refactorization is required online.



