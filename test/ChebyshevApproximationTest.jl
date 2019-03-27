using Plots
using Main.ChebyshevApproximation

function f(x)
    return exp.(-x.^2/0.1)
end

xFunc = -1:0.001:1
yFunc = f(xFunc)
plot(xFunc, yFunc)

@time Cns, C0 = GetCnList(f, 40)
a = [-1., -0.95, -0.96]
b = [0.1, 0.2, 0.3]
int = ChebyshevFitDefinedIntegral(a, b, Cns, C0)

fit = ChebyshevFit(xFunc, Cns, C0)

# fit = ChebyshevFitIntegral(0.123456789, Cns, C0) - ChebyshevFitIntegral(-0.123456789, Cns, C0)



plot!(xFunc, fit)
# plot(xFunc, fit .- yFunc .- 1)
