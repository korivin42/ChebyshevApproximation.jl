module ChebyshevApproximation
    export GetCnList, ChebyshevFit, ChebyshevFitIntegral, ChebyshevFitDefinedIntegral

    """
        Function for calculation Chebyshev polinomial of n-th order in a single point z or in array of points
        Tₙ(z) = cos( n ArcCos(z) )
    """
    function ChebyshevT(z, n::Int64)
        all(abs.(z) .<= 1) || error("Chebyshev T_n is defined on [-1, 1], value of |z|>1 was provided")
        return cos.(n .* acos.(z))
    end


    """
        Function for calculation integral of n-th Chebyshev polinomial in a single point z or in array of points
        Tₙ(z) = ∫ cos( n ArcCos(z) ) dz
    """
    function ChebyshevTintegral(z, n::Int64)
        all(abs.(z) .<= 1) || error("Chebyshev T_n is defined on [-1, 1], value of |z|>1 was provided")
        if n >= 2
            return n .* ChebyshevT(z, n+1) ./ (n^2 - 1) - z .* ChebyshevT(z, n) / (n - 1)
        elseif n ==1
            # return z.^2 / 2
            return ChebyshevT(z, 2) / 4 .+ 1//4
        end
    end


    """
        Function for calculation the n-th coefficient C_n in Chebyshev series
        f(z) = C₀/2 + Σₙ₌₁ᴺ Cₙ ⋅ Tₙ(z)

        Input: f -- function which is to be approximated,
               n -- index of the coefficient ([0..N])
               N -- maximal index of the series. ~20 is usually suffitient for the most of functions
    """
    function Cn(f::Function, n::Int64, N::Int64)
        ks = 1:1:N
        nodes = f(cos.(π * (ks .- 1/2) / N))
        coss = cos.(π * n * (ks .- 1/2) / N)
        return 2 / N * sum(nodes .* coss)
    end


    """
        Function to get a list of coefficients up to Nₘₐₓ-th order for Chebyshev approximation

        Input: f     -- function which is to be approximated,
               N_max -- maximal index of the series. ~20 is usually suffitient for the most of functions.
                        default is 20
        Output: Cns  -- List of coefficients for n ∈ [1, N_max]
                C0   -- 0-th order coeffitcent

    """
    function GetCnList(f::Function, Nmax::Int64=20)
        Cns = zeros(Nmax) .- NaN
        for n in 1:1:Nmax
            Cns[n] = Cn(f, n, Nmax)
        end
        C0 = Cn(f, 0, Nmax)
        return Cns, C0
    end


    """
        Function provides Chebyshev approximation for given coefficients
    """
    function ChebyshevFit(x, Cns, C0)
        func = zeros(length(x))
        for k in 1:1:length(Cns)
            func .+= Cns[k] .* ChebyshevApproximation.ChebyshevT(x, k)
        end
        func .+= C0 / 2.0
    end

    """
        Function provides integral of Chebyshev approximation for given coefficients
    """
    function ChebyshevFitIntegral(x, Cns, C0, integralC=0)
        func = zeros(length(x))
        for k in 1:1:length(Cns)
            func .+= Cns[k] .* ChebyshevTintegral(x, k)
        end
        func .+= C0 / 2.0 .* x .+ integralC
    end

    """
        Function provides defined integral ∫ₐᵇ of Chebyshev approximation for given coefficients
    """
    function ChebyshevFitDefinedIntegral(a, b, Cns, C0)
        ChebyshevFitIntegral(b, Cns, C0) .- ChebyshevFitIntegral(a, Cns, C0)
    end

end  # module ChebyshevApproximation
