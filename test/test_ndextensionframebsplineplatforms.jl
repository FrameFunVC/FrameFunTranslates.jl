
using FrameFunTranslates, Test, LinearAlgebra
@testset "approximation power" begin
    f = (x,y)->exp(x*y)
    for d1 in 1:2, d2 in 1:2, PLATFORM in (NdEpsBSplinePlatform, NdBSplinePlatform, NdCDBSplinePlatform), N1 in (10,15), N2 in (10,15)
        P = ExtensionFramePlatform(PLATFORM((d1,d2)), (0.0..0.5)^2); N = (N1,N2)
        plunge = plungeoperator(P,N); A = AZ_A(P,N); Zt = AZ_Zt(P,N)
        S = reducedAZ_AAZAreductionsolver(P,N;solverstyle=ReducedAZStyle())
        b = samplingoperator(P,N)*f
        x1 = S*plunge*b
        x2 = Zt*(b-A*x1)
        @test norm(A*(x1+x2)-b) < 5e-3
    end
end




using FrameFunTranslates, Test, LinearAlgebra
@testset "ExtensionFramePlatform, ReductionSolver AZ approximation power, Nd" begin
    f = (x,y)->exp(x*y)
    for d in 1:5, PLATFORM in (NdEpsBSplinePlatform, NdBSplinePlatform, NdCDBSplinePlatform)
        P = ExtensionFramePlatform(PLATFORM((d,d)), (0.0..0.5)^2); N = (30,30)
        F = Fun(f, P, N;solverstyle=ReducedAZStyle())
        b = samplingoperator(P,N)*f
        A = AZ_A(P,N)
        @test norm(A*coefficients(F)-b) < 5e-3
    end
end
