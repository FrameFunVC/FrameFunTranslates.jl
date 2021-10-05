
using Test
using LowRankApprox, LinearAlgebra, DomainSets, FrameFun

using FrameFunTranslates

@testset "AbstractBSplinePlatforms" begin
    include("test_abstractbsplineplatforms.jl")
end

@testset "ExtensionFrame BSpline Platform" begin
    include("test_extensionframebsplineplatforms.jl")
end

@testset "Nd BSpline Platforms" begin
    include("test_ndbsplineplatforms.jl")
end

@testset "Nd ExtensionFrame BSpline Platform" begin
    include("test_ndextensionframebsplineplatforms.jl")
end

@testset "nonzero_coefficients" begin
    include("test_nonzero_coefficients.jl")
end

@testset "1D FrameFunTranslates" begin
    include("test_1dbsplineextension.jl")
end

@testset "CDPET platform" begin
    include("test_compactperiodicequispacedtranslatesdual.jl")
end

@testset "AZ sparse" begin
    include("test_azsparse.jl")
end
