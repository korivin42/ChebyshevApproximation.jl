# ChebyshevApproximation.jl

A simple rootine for calculating 1D Chebyshev Approximation and its defined and undefined integrals

## General theory
A 1D Chebyshev Approximation means representing a function f(x) as a series of Chebyshev Polynomials of the first kind.
Note, that x should be from an interval [-1, 1]

## Examples:
Simple fitting:
```
  using ChebyshevApproximation          # loading the module
  f(x) = exp.(-x.^2/0.1)                # a random function defined x ∈ [-1, 1]
  Cns, C0 = GetCnList(f)                # calculating fitting coeffitients for the series

  xFit = -1:0.0001:1
  fit = ChebyshevFit(xFit, Cns, C0)     # obtaining fitted function in xFit points
```

Note, that `GetCnList` takes in a function `f` and optionally one can set the ammount of polynomials used in the approxiations as:
```
  GetCnList(f, 400)
```
here, we used 400 polynomials to fit function `f`. The default is 30.

Integration:
To obtain an exact value of defined integral of the fit (defined by a set of `Cns, C0`) one can run:
```
  ChebyshevFitDefinedIntegral(-1., 1., Cns, C0)
```
here `ChebyshevFitDefinedIntegral` takes in -1 as a lower limit and 1 as an upper one. Fit is defined by the set of `Cns, C0` coefficients.
