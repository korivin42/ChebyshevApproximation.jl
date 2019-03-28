# ChebyshevApproximation.jl

A simple rootine for calculating 1D Chebyshev Approximation and its defined and undefined integrals

## General theory
A 1D Chebyshev Approximation means representing a function ![equation](https://latex.codecogs.com/gif.latex?f%28z%29%20%5Capprox%20%5Cfrac%7Bc_0%7D%7B2%7D%20&plus;%20%5Csum_%7Bk%3D1%7D%5E%7BN%7D%20c_kT_k%28z%29) 

Example:
```
  f(x) = exp.(-x.^2/0.1)
```
