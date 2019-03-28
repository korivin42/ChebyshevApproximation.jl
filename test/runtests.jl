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

    @testset "cos" begin
        f(x) = cos.(x)
        Cns, C0 = GetCnList(f)
        for xFit in xTestList
            @test ChebyshevFit(xFit, Cns, C0)[1] ≈ f(xFit)
        end
    end

    @testset "cos(150x)^2" begin
        f(x) = cos.(150x).^2
        Cns, C0 = GetCnList(f, 400)
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


@testset "Testing integrating capapibilities" begin

    @testset "exp" begin
        f(x) = exp.(x)
        Cns, C0 = GetCnList(f)
        @test ChebyshevFitDefinedIntegral(-1., 1., Cns, C0)[1] ≈ f(1.) - f(-1.)
    end

    @testset "cos(150x)^2" begin
        f(x) = cos.(150x).^2
        F(x) = x./2 .+ sin.(300x) ./ 600
        Cns, C0 = GetCnList(f, 400)
        @test ChebyshevFitDefinedIntegral(-1., 1., Cns, C0)[1] ≈ F(1) - F(-1)
    end

    @testset "Gaussian" begin
        f(x) = exp.(-x.^2/0.1)
        Cns, C0 = GetCnList(f, 40)
        @test ChebyshevFitDefinedIntegral(-1., 1., Cns, C0)[1] ≈ 0.5604947810132863
    end

end
