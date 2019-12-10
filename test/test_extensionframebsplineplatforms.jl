
using FrameFunTranslates, Test, LinearAlgebra
@testset "approximation power" begin
    for d in 2:5
        for PLATFORM in (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform)
            P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); N = 30
            plunge = plungeoperator(P,N;L=4N); A = AZ_A(P,N;L=4N); Zt = AZ_Zt(P,N;L=4N)
            M = plunge*A; S = reducedAZ_AAZAreductionsolver(P,N;solverstyle=ReducedAZStyle())
            b = samplingoperator(P,N;L=4N)*exp
            x1 = S*plunge*b
            x2 = Zt*(b-A*x1)
            @test norm(A*(x1+x2) -b) < 8e-6
        end
    end
end

using FrameFunTranslates, Test
@testset "AZ approximation power" begin
    for d in 1:5
        for PLATFORM in (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform)
            P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); N = 30
            A = AZ_A(P,N)
            b = samplingoperator(P,N;L=4N)*exp
            F = Fun(exp, P, N;L=4N, solverstyle=ReducedAZStyle())
            @test norm(A*coefficients(F) -b) < 5e-3
        end
    end
end
