
using InfiniteVectors, CompactTranslatesDict.CompactInfiniteVectors
using GridArrays.ModCartesianIndicesBase: ModCartesianIndices
using BasisFunctions: Dictionary, Measure, AbstractGrid, grid, support
using GridArrays: mask


export nonzero_coefficient
"""
    nonzero_coefficient(dict::Dictionary, measure::Measure)

The coefficient indices of `dict` that overlap with the boundary of the measure.
"""
nonzero_coefficient(dict::Dictionary, measure::Measure)  =
    nonzero_cols(dict, measure)
nonzero_cols(dict::Dictionary, measure::Measure) =
    nonzero_cols(dict, grid(measure))
function nonzero_cols(dict::Dictionary, grid::AbstractGrid)
    q = div.(size(supergrid(grid)),size(dict))
    r = rem.(size(supergrid(grid)),size(dict))
    @assert all(r.==0)
    nonzero_cols(compactinfinitevector(dict, grid), q, mask(grid))
end

nonzero_cols(b::CompactInfiniteVector, m::Tuple{Int}, gridmask::BitArray{1}) =
    (s=InfiniteVectors.support(b);nonzero_cols(CartesianIndex(first(s)):CartesianIndex(last(s)), m, gridmask))

nonzero_cols(b::NTuple{N,CompactInfiniteVector}, m::NTuple{N,Int}, gridmask::BitArray{N}) where N =
    nonzero_cols(support(b), m, gridmask)

function nonzero_cols(b_support::CartesianIndices{N}, m::NTuple{N,Int}, gridmask::BitArray{N}) where N
    gridsize = size(gridmask)
    dictsize = div.(gridsize,m)
    A = falses(dictsize)
    for k in CartesianIndices(dictsize)
        l = (k.I .- 1).*m
        bk_support = b_support .+ CartesianIndex(l)
        support_in = false
        support_out = false
        for i in ModCartesianIndices(gridsize, first(bk_support), last(bk_support))
            if gridmask[i]
                support_in = true
            else
                support_out = true
            end
            if support_in && support_out
                A[k] = true
                break;
            end
        end
    end
    findall(A)
end
