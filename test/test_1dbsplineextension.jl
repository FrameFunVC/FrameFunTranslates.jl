
using FrameFunTranslates, Test, DomainSets, FrameFun
@testset "truncated size" begin
    Ns = 20:20:300
    ds = 1:4
    PLATFORMs = (BSplinePlatform, EpsBSplinePlatform, CDBSplinePlatform)
    crop_tols = 10.0.^(-16.:6.:-10.)
    colsizes = Array{Int}(undef, length(PLATFORMs), length(ds), length(Ns), length(crop_tols))
    rowsizes = similar(colsizes)
    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns), (l,crop_tol) in enumerate(crop_tols)
        P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5);
        colsizes[i,j,k,l], rowsizes[i,j,k,l] = size(reducedAAZAoperator(P,N;solverstyle=ReducedAZStyle(),nz_tol=crop_tol))
    end

    @test all(rowsizes[:,1,:,:] .== 2)
    @test all(rowsizes[:,2,:,:] .== 6)
    @test all(rowsizes[:,3,:,:] .== 6)
    @test all(rowsizes[:,4,:,:] .== 10)

    # Test CDBSplinePlatform
    @test all(0 .<= colsizes[3,1,:,1] .<= 7)
    @test all(0 .<= colsizes[3,1,:,2] .<= 0)
    @test all(19 .<= colsizes[3,2,:,1] .<= 24)
    @test all(12 .<= colsizes[3,2,:,2] .<= 12)
    @test all(20 .<= colsizes[3,4,:,2] .<= 20)
    @test all(33 .<= colsizes[3,4,:,1][2:end] .<= 42)

    @test all(176 .<= colsizes[2,1,:,1][end-1:end] .<= 176)
    @test all(122 .<= colsizes[2,1,:,2][end-1:end] .<= 122)
    @test all(316 .<= colsizes[2,2,:,1][end-1:end] .<= 316)
    @test all(196 .<= colsizes[2,2,:,2][end-1:end] .<= 196)
    @test all(308 .<= colsizes[2,4,:,2][end-1:end] .<= 308)
    @test all(544 .<= colsizes[2,4,:,1][2:end][end-1:end] .<= 546)

    @test all(122 .<= colsizes[1,1,:,2][end-4:end] .<= 122)
    @test all(196 .<= colsizes[1,2,:,2][end-4:end] .<= 196)
    @test all(254 .<= colsizes[1,3,:,2][end-4:end] .<= 254)
    @test all(308 .<= colsizes[1,4,:,2][end-4:end] .<= 308)
end

using FrameFunTranslates, Test
@testset "loss of information" begin
    Ns = [300,]
    ds = 1:4
    PLATFORMs = (BSplinePlatform, EpsBSplinePlatform, CDBSplinePlatform)
    crop_tols = [1e-8,]
    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns), (l,crop_tol) in enumerate(crop_tols)
        P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5);
        plunge = plungeoperator(P,N;L=4N); A = AZ_A(P,N;L=4N); Zt = AZ_Zt(P,N;L=4N);
        M = plunge*A;
        S = reducedAAZAoperator(P,N;solverstyle=ReducedAZStyle(),nz_tol=crop_tol)
        @test all(size(M) .> size(S))
        @test norm(M)≈norm(S)
    end
end


using Test, FrameFunTranslates, Statistics
@testset "errors, and timings" begin
    PLATFORMs = (EpsBSplinePlatform, BSplinePlatform, CDBSplinePlatform)
        Ns1 = [1<<k for k in 4:10]
        Ns2 = [1<<k for k in 9:16]
        ds = 1:4
        errors = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns1))
        timings = Array{Float64}(undef, length(PLATFORMs), length(ds), length(Ns2))

    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns1)
        P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5);
        F,_ = @timed Fun(exp, P, N;L=4N, solverstyle=ReducedAZStyle(),nz_tol=1e-10)
        errors[i,j,k] = abserror(exp, F)
    end

    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds), (k,N) in enumerate(Ns2)
        P = ExtensionFramePlatform(PLATFORM(d), 0.0..0.5);
        timings[i,j,k]= median([@timed(Fun(exp, P, N;L=4N, solverstyle=ReducedAZStyle(), nz_tol=1e-10))[2] for l in 1:4])
    end

    # Test if errors go down fast enough
    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds)
        @test all(errors[i,j,:] .<= max.(Float64.(Ns1).^(-d)*errors[i,j,1]*Ns1[1]^d, 1e-10))
    end

    # Test if method is fast enough
    for (i,PLATFORM) in enumerate(PLATFORMs), (j,d) in enumerate(ds)
        @test all(timings[i,j,:] .<= 5e-5Ns2)
    end

end
