
using FrameFunTranslates, Test
using FrameFunTranslates.CompactAZ.ReducedAZ.ReductionSolvers: ReductionSolver
@testset "approximation power" begin
    for d in 1:5
        for PLATFORM in (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform)
            P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); N = 30
            plunge = plungeoperator(P,N;L=4N); A = AZ_A(P,N;L=4N); Zt = AZ_Zt(P,N;L=4N)
            M = plunge*A; S = ReductionSolver(M;crop=true)
            b = samplingoperator(P,N;L=4N)*exp
            x1 = S*plunge*b
            x2 = Zt*(b-A*x1)
            F = DictFun(dictionary(P,N), x1 + x2)
            @test abserror(exp, F) < 5e-3
        end
    end
end

using FrameFunTranslates, Test
using FrameFunTranslates.CompactAZ.ReducedAZ.ReductionSolvers: ReductionSolver
@testset "AZ approximation power" begin
    for d in 1:5
        for PLATFORM in (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform)
            P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5); N = 30
            F = Fun(exp, P, N;REG=ReductionSolver, L=4N, crop=true)
            @test abserror(exp, F) < 5e-3
        end
    end
end
