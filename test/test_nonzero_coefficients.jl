
using FrameFunTranslates, Test
using FrameFunTranslates.CompactAZ.CompactFrameFunExtension: nonzero_rows, nonzero_coefficients
@testset "ExtensionFramePlatform, nonzero_rows, nonzero_coefficients" begin
    # 1D nonzero_coefficients
    for d in (1,2,3,4)
        P = ExtensionFramePlatform(EpsBSplinePlatform(d), 0.0..0.5)
        N = 100
        A = Matrix(AZ_A(P,N; L=4N))
        Zt = Matrix(AZ_Zt(P,N;L=4N))
        gt = nonzero_coefficients(P,N;L=4N)
        n = findall(sum(abs.(Matrix(A*Zt*A-A)); dims=1)[:] .> 1e-10)
        for ni in n
            @test ni ∈ gt
        end
    end

    # 1D nonzero_rows
    P = ExtensionFramePlatform(CDBSplinePlatform(), 0.0..0.5)
    N = 100
    M = plungeoperator(P,N;L=4N)*AZ_A(P,N;L=4N);size(M)
    @test findall(nonzero_rows(Matrix(M),nz_tol=1e-10))[:] ==
        findall(sum(abs.(Matrix(M));dims=2)[:] .> 1e-10)

    # 2D nonzero_coefficients
    P = ExtensionFramePlatform(NdEpsBSplinePlatform((3,3)),(0.0..0.5)^2)
    N = (10,10)
    A = Matrix(AZ_A(P,N))
    Zt = Matrix(AZ_Zt(P,N))
    gt = nonzero_coefficients(P,N)
    n = findall(reshape(sum(abs.(Matrix(A*Zt*A-A)); dims=1),N) .> 1e-10)
    for ni in n
        @test ni ∈ gt
    end

    # 2D nonzero_rows
    P = ExtensionFramePlatform(NdCDBSplinePlatform((3,3)),(0.0..0.5)^2)
    N = (10,10)
    M = plungeoperator(P,N)*AZ_A(P,N);size(M)
    @test findall(nonzero_rows(Matrix(M),nz_tol=1e-10))[:] ==
        findall(sum(abs.(Matrix(M));dims=2)[:] .> 1e-10)
end
