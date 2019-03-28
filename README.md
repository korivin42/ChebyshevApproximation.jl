# ChebyshevApproximation.jl

A simple rootine for calculating 1D Chebyshev Approximation and its defined and undefined integrals

## General theory
A 1D Chebyshev Approximation means representing a function f(x) as a series of Chebyshev Polynomials of the first kind.
Note, that x should be from an interval [-1, 1]

Example:
```
  using ChebyshevApproximation          # loading the module
  f(x) = exp.(-x.^2/0.1)                # a random function defined x ∈ [-1, 1]
  Cns, C0 = GetCnList(f)                # calculating fitting coeffitients for the series

  xFit = -1:0.0001:1
  fit = ChebyshevFit(xFit, Cns, C0)     # obtaining fitted function in xFit points
```
