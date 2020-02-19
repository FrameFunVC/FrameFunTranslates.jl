
using LinearAlgebra, FrameFunTranslates, Test, DomainSets, FrameFun
using FrameFun.ApproximationProblems: approximationproblem

@testset "(dual)dictionaries" begin
    P = BSplinePlatform()
    B = dictionary(P,10)

    ap = approximationproblem(P,10)
    μ1 = discretemeasure(SamplingStyle(ap), ap)
    μ2 = discretemeasure(sampling_grid(P,10))
    @test μ1≈μ2


    d1 = dictionary(P,10)
    d2 = azdual_dict(P,10)

    g1 = mixedgram(d1, d2)
    g2 = mixedgram(d1, d2, discretemeasure(sampling_grid(P,10)))
    @test g1 isa CirculantOperator
    @test g2 isa CirculantOperator
    @test g2 ≈ IdentityOperator(B)

    P = EpsBSplinePlatform()
    g = sampling_grid(P,10)
    d1 = dictionary(P,1000)
    d2 = azdual_dict(P,1000;threshold=1e-6)
    @test (operator(d2) isa VerticalBandedOperator)
    g2 = mixedgram(d1, d2, discretemeasure(sampling_grid(P,1000)))
    @test norm(IdentityOperator(d1)-g2) < 1e-4


    P = CDBSplinePlatform(5)
    d1 = dictionary(P,20)
    d2 = azdual_dict(P,20)
    @test d2 isa FrameFunTranslates.BSplinePlatforms.CompactPeriodicEquispacedTranslatesDuals.CompactPeriodicEquispacedTranslatesDual
    g2 = mixedgram(d1, d2, discretemeasure(sampling_grid(P,20)))
    @test IdentityOperator(d1)≈g2


    P = EpsBSplinePlatform()
    @test 20==length(sampling_grid(P,10; oversamplingfactor=1.6,samplingstyle=OversamplingStyle()))
    @test 20==FrameFun.samplingparameter(P,10; oversamplingfactor=1.6,samplingstyle=OversamplingStyle())
end

using Test, FrameFunTranslates
@testset "samplng parameter, oversampling" begin
    P = BSplinePlatform()
    for N in 1:100
        @test divrem(samplingparameter(P, N; samplingstyle=OversamplingStyle(), oversamplingfactor=pi) ,N)[2] == 0
    end

    P = ExtensionFramePlatform(BSplinePlatform(), 0.0..0.5)
    for N in 1:100
        @test divrem(samplingparameter(P, N; samplingstyle=OversamplingStyle(), oversamplingfactor=pi) ,N)[2] == 0
    end

    P = NdBSplinePlatform((3,3))
    for N in 1:100
        @test (0,0) == map(x->divrem(x, N)[2], samplingparameter(P, (N,N); samplingstyle=ProductSamplingStyle(OversamplingStyle(),OversamplingStyle()), oversamplingfactor=pi))
    end

    P = ExtensionFramePlatform(NdBSplinePlatform((3,3)), (0.0..0.5)^2)
    for N in 5:100
        @test (0,0) == map(x->divrem(x, N)[2], samplingparameter(P, (N,N); oversamplingfactor=pi))
    end


    P = ExtensionFramePlatform(NdBSplinePlatform((1,3)),(0.0..0.5)^2)
    N = 10
    @test all(operator.(elements(basis(azdual_dict(P,(N,N);L=(4N,4N))))) .≈
        operator.(elements(basis(azdual_dict(P,(N,N),oversamplingfactor=4)))))
end

using Test, FrameFunTranslates
@testset "Basis platform approximation" begin

    P1 = BSplinePlatform()
    P2 = EpsBSplinePlatform()
    P3 = CDBSplinePlatform()
    @test SamplingStyle(approximationproblem(P1,10)) == InterpolationStyle()
    @test SamplingStyle(approximationproblem(P2,10)) == InterpolationStyle()
    @test SamplingStyle(approximationproblem(P3,10)) == OversamplingStyle()

    # Test the approximation power
    f = x->exp(cos(10pi*x))
    N = 60
    x = PeriodicEquispacedGrid(10N, support(dictionary(P1,10)))
    t = .123
    for P in (P1,P2,P3)
        @test SolverStyle(InterpolationStyle(), approximationproblem(P, 10)) == DualStyle()
        F = Fun(f, P, N)
        @test norm(F.(x)-f.(x),Inf)  < .007
        @test abs(F(t)-f(t)) < .0005
    end
end


using FrameFunTranslates, Test, LowRankApprox
@testset "ExtensionFramePlatform, ReductionSolver approximation power" begin
    for d in 1:5
        for PLATFORM in (BSplinePlatform, EpsBSplinePlatform, CDBSplinePlatform)
            P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); N = 30
            plunge = plungeoperator(P,N;L=4N); A = AZ_A(P,N;L=4N); Zt = AZ_Zt(P,N;L=4N)
            S = reducedAZ_AAZAreductionsolver(P,N;solverstyle=ReducedAZStyle(),lraoptions=LRAOptions(atol=1e-14))
            b = samplingoperator(P,N;L=4N)*exp
            x1 = S*plunge*b
            x2 = Zt*(b-A*x1)
            @test norm(A*(x1+x2)-b) < 5e-3
        end
    end
end
