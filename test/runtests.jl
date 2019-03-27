using ChebyshevApproximation
using Test

@testset "Testing fitting/interpolation capapibilities" begin
    xTestList = -1.:0.2:1.

    @testset "exp" begin
        f(x) = exp.(x)
        Cns, C0 = GetCnList(f)
        for xFit in xTestList
            @test ChebyshevFit(xFit, Cns, C0)[1] ≈ f(xFit)
        end
    end

    @testset "sin" begin
        f(x) = cos.(x)
        Cns, C0 = GetCnList(f)
        for xFit in xTestList
            @test ChebyshevFit(xFit, Cns, C0)[1] ≈ f(xFit)
        end
    end

    @testset "Gaussian" begin
        f(x) = exp.(-x.^2/0.1)
        Cns, C0 = GetCnList(f, 40)
        for xFit in xTestList
            @test ChebyshevFit(xFit, Cns, C0)[1] ≈ f(xFit)
        end
    end
end
