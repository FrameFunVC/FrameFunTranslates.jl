

using Test, BasisFunctions, FrameFun, FrameFunTranslates, CompactTranslatesDict, CardinalBSplines
using GridArrays: similargrid

N = 10

x = LinRange(-3,3,100)

B = GenericPeriodicEquispacedTranslates(PeriodicEquispacedGrid(N,-1.123,1.432), x->CenteredBSpline(2)(N*x),(-2..2)/N)

m = 2
D = CompactPeriodicEquispacedTranslatesDual(B, m)
m1 = evaluation(B,similargrid(interpolation_grid(B), Float64,m*N))
m2 = evaluation(D,similargrid(interpolation_grid(B), Float64,m*N))
@test m1'm2≈IdentityOperator(B, B)
P = platform(B)
@test P isa CDPETPlatform
μ = discretemeasure(similargrid(interpolation_grid(B),Float64,m*N))
D = dualdictionary(P, N, μ)


@test mixedgram(B,D,μ)≈IdentityOperator(B, B)




P = CDPETPlatform(BSplineTranslatesBasis(6,3,-1,1))
d1 = dictionary(P,6)
d2 = azdual(P,6)
g2 = mixedgram(d1, d2, discretemeasure(sampling_grid(P,6)))

using InfiniteVectors
b = CompactTranslatesDict.signal(P.dict,2)
primal_signal = PeriodicInfiniteVector(b, 12)[0:11]
c = inv(b, 2,K=d2.minimalK)
dual_signal = PeriodicInfiniteVector(c, 12)[0:11]

@test evaluation(d1, sampling_grid(P,6)).A[:,1]≈primal_signal
@test evaluation(d2, sampling_grid(P,6)).A[:,1]≈dual_signal

@test g2≈IdentityOperator(d1)
