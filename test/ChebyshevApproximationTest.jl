using Plots
using Main.ChebyshevApproximation

function f(x)
    return exp.(x)
end

xFunc = -1:0.001:1
yFunc = f(xFunc)
plot(xFunc, yFunc)

@time Cns, C0 = GetCnList(f)
a = 0.123456789
b = 0.56789101112
int = ChebyshevApproximation.ChebyshevFitDefinedIntegral(a, b, Cns, C0)

# fit = ChebyshevFit(xFunc, Cns, C0)

# fit = ChebyshevFitIntegral(0.123456789, Cns, C0) - ChebyshevFitIntegral(-0.123456789, Cns, C0)



# plot!(xFunc, fit)
# plot(xFunc, fit .- yFunc .- 1)
