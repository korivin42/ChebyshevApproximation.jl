# ChebyshevApproximation.jl

A simple rootine for calculating 1D Chebyshev Approximation and its defined and undefined integrals

## General theory
A 1D Chebyshev Approximation means representing a function ![equation](https://latex.codecogs.com/gif.latex?f(z)&space;\approx&space;\frac{c_0}{2}&space;&plus;&space;\sum_{k=1}^{N}&space;c_kT_k(z)) 

Example:
```
  f(x) = exp.(-x.^2/0.1)
```
