module BSplinePlatforms

using FrameFun.Platforms, FrameFun.BasisFunctions, CompactTranslatesDict
using FrameFun.BasisFunctions: AbstractMeasure

import FrameFun.Platforms: dictionary, SolverStyle, measure, SamplingStyle,
    dualdictionary, correctparamformat, unsafe_dictionary, platform
import FrameFun.FrameFunInterface: correct_sampling_parameter, regularization_threshold, SamplingStrategy


abstract type AbstractPeriodicEquispacedTranslatesPlatform{T,S} <: BasisPlatform end

SolverStyle(p::AbstractPeriodicEquispacedTranslatesPlatform, ::SamplingStyle) = DualStyle()
correct_sampling_parameter(::SamplingStrategy, ::AbstractPeriodicEquispacedTranslatesPlatform, param, L; options...) = error()
correct_sampling_parameter(::SamplingStrategy, ::AbstractPeriodicEquispacedTranslatesPlatform, param::Int, L::Int; options...) =
    (round(Int, L/param) * param)
correctparamformat(::AbstractPeriodicEquispacedTranslatesPlatform, ::Int) = true

abstract type AbstractEpsPeriodicEquispacedTranslatesPlatform{T,S} <: AbstractPeriodicEquispacedTranslatesPlatform{T,S} end

function dualdictionary(platform::AbstractEpsPeriodicEquispacedTranslatesPlatform{T}, param, measure::AbstractMeasure;
        threshold=regularization_threshold(T), options...) where T
    dual = gramdual(dictionary(platform, param), measure;options...)
    op = BasisFunctions.operator(dual)
    @assert op isa CirculantOperator

    # Get the abs of the first row to determine the bandwidth
    e = zeros(src(op)); e[1] = 1
    column = op*e
    # bandwidth determined by threshold
    # we want \|I-G (G̃+E)\|< threshold, therefore \|E\| < threshold / \|G\|
    mask = abs.(column) .< threshold / norm(inv(op).A.vcvr_dft,Inf)
    bw = (findfirst(mask), findlast(mask))
    if bw[1] == nothing
        return dual
    end
    bw = bw[1]-1,bw[2]+1
    a = [column[bw[2]:end];column[1:bw[1]]]
    # Create bandlimited operator (which is also circulant)
    op_replace = VerticalBandedOperator(src(op), dest(op), a, 1, -1+bw[2]-length(dual))
    op_replace*src(op)
end

abstract type AbstractCDPeriodicEquispacedTranslatesPlatform{T,S}<: AbstractPeriodicEquispacedTranslatesPlatform{T,S} end
SamplingStyle(::AbstractCDPeriodicEquispacedTranslatesPlatform) = OversamplingStyle()


dualdictionary(platform::AbstractCDPeriodicEquispacedTranslatesPlatform, param, measure::AbstractMeasure; options...) =
    error("No azdual_dict for `CDBSplinePlatform` and $(typeof(measure))")

function dualdictionary(platform::AbstractCDPeriodicEquispacedTranslatesPlatform, param, measure::UniformDiracComb;
        options...)
    dict = dictionary(platform, param)
    g = points(measure)
    @assert support(dict) ≈ coverdomain(g)
    m = length(g) / length(dict)
    @assert round(Int,m) ≈ m
    m = round(Int, m)
    if m == 1
        @warn "No compact dual possible, try oversampling"
        return dual(dict, measure; options...)
    else
        if isperiodic(g)
            return CompactPeriodicEquispacedTranslatesDual(dict, m; options...)
        else
            error()
        end
    end
end


export BSplinePlatform
"""
    struct BSplinePlatform{T,D} <: AbstractPeriodicEquispacedTranslatesPlatform{T,T}

A platform of equispaced periodic translates of B-spline of a given order.
The dual dictionaries are determined by the inverse of the Gram matrix.


See also: [`EpsBSplinePlatform`](@ref), [`CDBSplinePlatform`](@ref)

# Example
```jldocs
julia> P = BSplinePlatform()
BSplinePlatform{Float64,3}()

julia> g = sampling_grid(P,10)
10-element PeriodicEquispacedGrid{Float64}:
 0.0
 0.09999999999999999
 0.19999999999999998
 0.3
 ⋮
 0.7000000000000001
 0.7999999999999999
 0.9

julia> d1 = dictionary(P,10)
Periodic equispaced translates of B spline of degree 3
    ↳ length = 10
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)



julia> d2 = azdual_dict(P,10)
Dictionary M * P

P   :   Periodic equispaced translates of B spline of degree 3
            ↳ length = 10
            ↳ Float64 -> Float64
            ↳ support = 0.0..1.0 (Unit)
M   :   Multiplication by Circulant{Float64,Complex{Float64}}


julia> mixedgram(d1,d2,discretemeasure(g))≈IdentityOperator(d1,d2)
true
```
"""
struct BSplinePlatform{T,D} <: AbstractPeriodicEquispacedTranslatesPlatform{T,T}
end

BSplinePlatform(degree::Int=3) = BSplinePlatform{Float64,degree}()
unsafe_dictionary(platform::BSplinePlatform{T,D}, param::Int) where {T,D} = BSplineTranslatesBasis{T,D,false}(param)


export EpsBSplinePlatform
"""
    struct EpsBSplinePlatform{T,D} <: AbstractEpsPeriodicEquispacedTranslatesPlatform{T,T}

A platform of equispaced periodic translates of B-spline of a given order.
The dual dictionaries are determined by the inverse of the Gram matrix. This Gram
matrix has exponentially decaying elements. This platform truncates elemens lower
`threshold` which results in a banded Gram.

See also: [`BSplinePlatform`](@ref), [`CDBSplinePlatform`](@ref)

# Example
```jldocs
julia> P = EpsBSplinePlatform()
EpsBSplinePlatform{Float64,3}()

julia> d1 = dictionary(P,1000)
Periodic equispaced translates of B spline of degree 3
    ↳ length = 1000
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)



julia> d2 = azdual_dict(P,1000;threshold=1e-4)
Dictionary M * P

P   :   Periodic equispaced translates of B spline of degree 3
            ↳ length = 1000
            ↳ Float64 -> Float64
            ↳ support = 0.0..1.0 (Unit)
M   :   Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}


julia> g2 = mixedgram(d1, d2, discretemeasure(sampling_grid(P,1000)))
Operator M₂ * M₁

M₂  :   Multiplication by Circulant{Float64,Complex{Float64}}
M₁  :   Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}


julia> norm(IdentityOperator(d1)-g2)
0.0004146509264381403
```
"""
struct EpsBSplinePlatform{T,D} <: AbstractEpsPeriodicEquispacedTranslatesPlatform{T,T}
end

EpsBSplinePlatform(degree::Int=3) = EpsBSplinePlatform{Float64,degree}()
unsafe_dictionary(platform::EpsBSplinePlatform{T,D}, param::Int) where {T,D} = BSplineTranslatesBasis{T,D,false}(param)

export CDBSplinePlatform
"""
    struct CDBSplinePlatform{T,D} <: AbstractCDPeriodicEquispacedTranslatesPlatform{T,T}

A platform of equispaced periodic translates of B-spline of a given order.
Their duals dictionaries are compact, but discrete, i.e.,
their values are known in an `PeriodicEquispacedGrid`

See also: [`DiscreteBSplineDict`](@ref), [`BSplinePlatform`](@ref), [`EpsBSplinePlatform`](@ref)

# Example
```jldocs
julia> P = CDBSplinePlatform(5)
CDBSplinePlatform{Float64,5}()

julia> d1 = dictionary(P,20)
Periodic equispaced translates of B spline of degree 5
    ↳ length = 20
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)
    ↳ degree = 5



julia> d2 = azdual_dict(P,20)
Equispaced translates of a discrete kernel dual to B-spline
    ↳ length = 20
    ↳ Float64 -> Float64
    ↳ support = 0.0..1.0 (Unit)
    ↳ degree = 5
    ↳ m = 2



julia> g2 = mixedgram(d1, d2, discretemeasure(sampling_grid(P,20)))
Operator M₂ * M₁

M₂  :   Multiplication by BasisFunctions.HorizontalBandedMatrix{Float64}
M₁  :   Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}


julia> IdentityOperator(d1)≈g2
true
```
"""
struct CDBSplinePlatform{T,D} <: AbstractCDPeriodicEquispacedTranslatesPlatform{T,T}
end

CDBSplinePlatform(degree::Int=3) = CDBSplinePlatform{Float64,degree}()
unsafe_dictionary(platform::CDBSplinePlatform{T,D}, param::Int) where {T,D} = BSplineTranslatesBasis{T,D,false}(param)

export GaussSplinePlatform, CDGaussSplinePlatform
struct GaussSplinePlatform{T,D} <: AbstractPeriodicEquispacedTranslatesPlatform{T,T}
end
GaussSplinePlatform(degree::Int=3) = GaussSplinePlatform{Float64,degree}()
unsafe_dictionary(platform::GaussSplinePlatform{T,D}, param::Int) where {T,D} = GaussTranslatesBasis{T,D,false}(param)

struct CDGaussSplinePlatform{T,D} <: AbstractCDPeriodicEquispacedTranslatesPlatform{T,T}
end
CDGaussSplinePlatform(degree::Int=3) = CDGaussSplinePlatform{Float64,degree}()
unsafe_dictionary(platform::CDGaussSplinePlatform{T,D}, param::Int) where {T,D} = GaussTranslatesBasis{T,D,false}(param)



using CompactTranslatesDict: PeriodicEquispacedTranslates
abstract type AbstractPETPlatorm{T,S,DICT} <: AbstractPeriodicEquispacedTranslatesPlatform{T,S} end
abstract type AbstractCDPETPlatorm{T,S,DICT} <: AbstractCDPeriodicEquispacedTranslatesPlatform{T,S} end
dictionary(platform::AbstractPETPlatorm, param) = similar(platform.dict, param)
dictionary(platform::AbstractCDPETPlatorm, param) = similar(platform.dict, param)
SamplingStyle(::AbstractPETPlatorm) = OversamplingStyle()
SamplingStyle(::AbstractCDPETPlatorm) = OversamplingStyle()

export PETPlatform
"""
    struct PETPlatform{T,S,DICT<:PeriodicEquispacedTranslates{T,S}} <: AbstractPETPlatorm{T,S,DICT}

A platform of periodic equispaced translates of a kernel.

See also: [`CDPETPlatform`](@ref), [`BSplinePlatform`](@ref)
# Example
```jldocs
julia> P = PETPlatform(BSplineTranslatesBasis(6,3,-1,1))
PETPlatform{Float64,Float64,GenericPeriodicEquispacedTranslates{Float64,Float64}}(Periodic equispaced translates of a periodic kernel function
    ↳ length = 6
    ↳ Float64 -> Float64
    ↳ support = -1.0..1.0

)

julia> d1 = dictionary(P,6)
Periodic equispaced translates of a periodic kernel function
    ↳ length = 6
    ↳ Float64 -> Float64
    ↳ support = -1.0..1.0



julia> d2 = azdual_dict(P,6)
Dictionary M * P

P   :   Periodic equispaced translates of a periodic kernel function
            ↳ length = 6
            ↳ Float64 -> Float64
            ↳ support = -1.0..1.0
M   :   Multiplication by Circulant{Float64,Complex{Float64}}


julia> g2 = mixedgram(d1, d2, discretemeasure(sampling_grid(P,6)))
Multiplication by Circulant{Float64,Complex{Float64}}



julia> g2≈IdentityOperator(d1)
true
```
"""
struct PETPlatform{T,S,DICT<:PeriodicEquispacedTranslates{T,S}} <: AbstractPETPlatorm{T,S,DICT}
    dict    :: DICT
end

export CDPETPlatform
"""
    struct CDPETPlatform{T,S,DICT<:PeriodicEquispacedTranslates{T,S}} <: AbstractCDPETPlatorm{T,S,DICT}

A platform of periodic equispaced translates of a kernel.
Their duals dictionaries are compact, but discrete, i.e.,
their values are known in an `PeriodicEquispacedGrid`. For some kernels these compact duals do not exist.

See also: [`PETPlatform`](@ref), [`CDBSplinePlatform`](@ref)
# Example
```jldocs
julia> P = CDPETPlatform(BSplineTranslatesBasis(6,3,-1,1))
CDPETPlatform{Float64,Float64,GenericPeriodicEquispacedTranslates{Float64,Float64}}(Periodic equispaced translates of a periodic kernel function
    ↳ length = 6
    ↳ Float64 -> Float64
    ↳ support = -1.0..1.0

)

julia> d1 = dictionary(P,6)
Periodic equispaced translates of a periodic kernel function
    ↳ length = 6
    ↳ Float64 -> Float64
    ↳ support = -1.0..1.0



julia> d2 = azdual_dict(P,6)
Equispaced translates of a discrete kernel dual
    ↳ Periodic equispaced translates of a periodic kernel function
      ↳ length = 6
      ↳ Float64 -> Float64
      ↳ support = -1.0..1.0
    ↳ m = 2



julia> g2 = mixedgram(d1, d2, discretemeasure(sampling_grid(P,6)))
Operator M₁ * M₂

M₂  :   Multiplication by BasisFunctions.VerticalBandedMatrix{Float64}
M₁  :   Multiplication by BasisFunctions.HorizontalBandedMatrix{Float64}


julia> g2≈IdentityOperator(d1)
true

```
"""
struct CDPETPlatform{T,S,DICT<:PeriodicEquispacedTranslates{T,S}} <: AbstractCDPETPlatorm{T,S,DICT}
    dict    :: DICT
end

# include("CompactPeriodicEquispacedTranslatesDuals.jl")
using .CompactPeriodicEquispacedTranslatesDuals


platform(dict::BSplineTranslatesBasis) =
    CDBSplinePlatform{coefficienttype(dict),degree(dict)}()
platform(dict::PeriodicEquispacedTranslates) =
    CDPETPlatform(dict)

end
