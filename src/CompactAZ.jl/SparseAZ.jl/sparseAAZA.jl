
using InfiniteVectors, GridArrays, FrameFun.BasisFunctions,
    CompactTranslatesDict, SparseArrays, FrameFun, CompactTranslatesDict.CompactInfiniteVectors

using ..ReducedAZ.ReductionSolvers: nonzero_cols, nonzero_rows
using FrameFun.BasisFunctions: DiscreteMeasure, support
using GridArrays.ModCartesianIndicesBase: ModCartesianIndices

function sparseAAZAmatrix(dict1::Dictionary, dict2::Dictionary, grid::AbstractGrid; nztol=0, opts...)
    ix1 = nonzero_cols(dict1, grid)
    A = sparseRAE(dict1, grid, ix1; opts...)
    # ix2 = nonzero_cols(dict2, grid)
    ix2 = overlappingindices(dict2, grid, findall(nonzero_rows(A)))
    ix3 = unique(sort(vcat(ix1,ix2)))
    Z = sparseRAE(dict2, grid, ix3; opts...)
    ImZA = sparseidentity(ix1,ix3)-Z'A
    droptol!(ImZA, nztol)
    RAE = sparseRAE(dict1, grid, ix3; opts...)
    droptol!(RAE*ImZA, nztol)
end

function overlappingindices(dict::Dictionary, grid::AbstractGrid, ix)
    boundarygrid_subindices = subindices(grid)[ix]
    if first(boundarygrid_subindices) isa Int
        boundarygrid_subindices = map(CartesianIndex, boundarygrid_subindices)
    end
    dict2mask = falses(size(dict))
    dict2support = support(compactinfinitevector(dict,grid))
    m = div.(size(supergrid(grid)),size(dict))
    for ix in boundarygrid_subindices
        coefindices =  CartesianIndex(cld.((ix-last(dict2support)).I,m).+1):CartesianIndex(fld.((ix-first(dict2support)).I,m).+1)
        if length(coefindices) > 0
            for l in ModCartesianIndices(size(dict2mask),first(coefindices),last(coefindices))
                dict2mask[l] = true
            end
        end
    end
    ix2 = findall(dict2mask)
end

function sparseidentity(ix1::Union{Vector{CartesianIndex{N}},Vector{Int}},ix2::Union{Vector{CartesianIndex{N}},Vector{Int}}) where N
    if length(ix1)==0 || length(ix2)==0
        return SparseMatrixCSC(length(ix2), length(ix1), ones(Int,length(ix1)+1), Int[], Int[])
    end
    R = Vector{Int}(undef, min(length(ix1),length(ix2)))
    colptr = Vector{Int}(undef, length(ix1)+1)
    colptr[1] = 1
    ix1i = 1
    ix2i = 1
    ix1e = first(ix1)
    ix2e = first(ix2)
    Ri = 0
    for (ix1i, ix1e) in enumerate(ix1)
        while ix2e < ix1e
            ix2i += 1
            if ix2i > length(ix2)
                break
            end
            ix2e=ix2[ix2i]
        end
        if ix1e == ix2e
            colptr[ix1i+1] = colptr[ix1i]+1
            Ri += 1
            R[Ri] = ix2i
        else
            colptr[ix1i+1] = colptr[ix1i]
        end
    end
    resize!(R,Ri)
    SparseMatrixCSC(length(ix2), length(ix1), colptr, R, ones(Int, Ri))
end


sparseRAE(b::CompactInfiniteVector,
    grid::AbstractGrid, indices::AbstractVector{Int},
    dictsize::NTuple{N,Int}; options...) where N =
    sparseRAE(tuple(b), grid, map(CartesianIndex,indices), dictsize; options...)

function sparseRAE(b::NTuple{N,CompactInfiniteVector},
    grid::AbstractGrid, indices::AbstractVector{CartesianIndex{N}},
    dictsize::NTuple{N,Int}; atol=1e-14) where N
    gridsize = size(supergrid(grid))

    B = Array(OuterProductArray(map(subvector, b)...))[:]
    B[abs.(B).<atol].=0
    b_support = support(b)
    b_supportlength = length(b_support)

    m = div.(gridsize,dictsize)
    L = LinearIndices(gridsize)

    nnz = b_supportlength*length(indices)

    colptr = Vector{Int}(undef, length(indices)+1)
    nzvals = Vector{eltype(B)}(undef, nnz)
    rowvals = Vector{Int}(undef, nnz)

    rowvalscol = Vector{Int}(undef, b_supportlength)
    nzvalscol = Vector{eltype(B)}(undef, b_supportlength)
    colix = Vector{Int}(undef, b_supportlength)

    gridmask = mask(grid)
    newindices = cumsum(gridmask[:])

    colptr[1] = 1
    nzvalindex = 1
    for (i,k) in enumerate(indices)
        # support of element with index k
        bk_support = CartesianIndex(m.*(k.I.-1)) .+ b_support
        bk_support_indices = ModCartesianIndices(gridsize, first(bk_support), last(bk_support))
        colptrcol = 0
        for (j,l) in enumerate(bk_support_indices)
            if gridmask[l] && B[j] != 0
                colptrcol += 1
                rowvalscol[colptrcol] = newindices[L[l]]
                nzvalscol[colptrcol] = B[j]
            end
        end
        for j in colptrcol+1:b_supportlength
            rowvalscol[j] = length(L)+1
        end

        for j in 1:b_supportlength
            colix[j] = j
        end
        sort!(colix ,1,b_supportlength, InsertionSort,Base.Perm(Base.Order.ForwardOrdering(),rowvalscol))
        # sortperm!(colix, rowvalscol)
        for j in 1:colptrcol
            rowvals[nzvalindex] = rowvalscol[colix[j]]
            nzvals[nzvalindex] = nzvalscol[colix[j]]
            nzvalindex += 1
        end
        colptr[i+1] = colptr[i] + colptrcol
    end

    SparseMatrixCSC(length(grid),length(indices),colptr,resize!(rowvals,nzvalindex-1),resize!(nzvals,nzvalindex-1))
end

sparseRAE(dict::Dictionary, grid::AbstractGrid, indices;opts...) =
    sparseRAE(compactinfinitevector(dict, grid), grid, indices, size(dict);opts...)



# sparsemixedgramcomplement(dict1::Dictionary, dict2::Dictionary, measure::DiscreteMeasure; opts...) =
#     sparsemixedgramcomplement(nonzero_cols(dict2, measure), compactinfinitevector(dict1, grid(measure)),
#         compactinfinitevector(dict2, grid(measure)), _oversamplingfactor(dict1,grid(measure)), grid(measure);opts...)
#
# _oversamplingfactor(dict::ExtensionFrame, grid::AbstractGrid) =
#     _oversamplingfactor(basis(dict), supergrid(grid))
# _oversamplingfactor(dict::TensorProductDict, grid::ProductGrid) =
#     map(_oversamplingfactor, elements(dict), elements(grid))
# function _oversamplingfactor(dict::Dictionary, grid::AbstractIntervalGrid)
#     q, r = divrem(length(grid), length(dict))
#     @assert r == 0
#     q
# end
#
# function difference_indices(b::NTuple{N,CompactInfiniteVector}, b̃::NTuple{N,CompactInfiniteVector}, m::NTuple{N,Int}) where N
#     b_support = support(b)
#     b̃_support = support(b̃)
#
#     dff1 = first(b_support)-last(b̃_support)
#     dff2 = last(b_support)-first(b̃_support)
#
#     myfun = (a,b) -> sign(a)*fld(abs(a), b)
#     CartesianIndex(map(myfun, dff1.I, m)):CartesianIndex(map(myfun, dff2.I,m))
# end
#
#
# function sparsemixedgramcomplement_nzband(indices::AbstractArray{CartesianIndex{N}}, T,
#         b̃::NTuple{N,CompactInfiniteVector}, b̃_support::CartesianIndices{N},
#         b::NTuple{N,CompactInfiniteVector}, b_support::CartesianIndices{N},
#         m::NTuple{N,Int}, dff::CartesianIndices{N}, grid::AbstractGrid) where N
#     L = length(dff)
#
#     R = Matrix{T}(undef, L, length(indices))
#     zerocartesian = CartesianIndex(ntuple(k->0,Val(N)))
#
#     bb = Array(OuterProductArray(map(subvector, b)...))
#     b̃b̃ = Array(OuterProductArray(map(subvector, b̃)...))
#
#     gridmask = mask(grid)
#     gridsize = size(supergrid(grid))
#
#     for (i,k) in enumerate(indices)
#         for (j,d) in enumerate(dff)
#             R[j,i] = - _innerproduct(bb, b̃b̃, b_support, b̃_support, m, k, d, gridmask, gridsize)
#             if d==zerocartesian
#                 R[j,i] += 1
#             end
#         end
#     end
#     R
# end
#
#
#
# function sparsemixedgramcomplement_nzband2(indices::AbstractArray{CartesianIndex{N}}, T,
#         b̃::NTuple{N,CompactInfiniteVector}, b̃_support::CartesianIndices{N},
#         b::NTuple{N,CompactInfiniteVector}, b_support::CartesianIndices{N},
#         m::NTuple{N,Int}, dff::CartesianIndices{N}, grid::AbstractGrid) where N
#     L = length(dff)
#
#     R = Matrix{T}(undef, L, length(indices))
#     zerocartesian = CartesianIndex(ntuple(k->0,Val(N)))
#
#     bb = Array(OuterProductArray(map(subvector, b)...))
#     b̃b̃ = Array(OuterProductArray(map(subvector, b̃)...))
#
#     gridmask = mask(grid)
#     gridsize = size(supergrid(grid))
#
#     L1 = length(bb)
#     nzvals1 = Vector{T}(undef, L1)
#     colptr1 = Vector{Int}(undef, 2)
#     rowvals1 = Vector{Int}(undef, L1)
#     v1 = SparseMatrixCSC(length(gridmask), 1, colptr1, rowvals1, nzvals1)
#
#     gridmaskcache1 = falses(size(bb))
#     rowvalscache1 = Vector{Int}(undef,L1)
#     nzvalscache1 = Vector{T}(undef,L1)
#     indexcache1 = Vector{Int}(undef,L1)
#
#     L2 = length(b̃b̃)
#     nzvals2 = Vector{T}(undef, L2)
#     colptr2 = Vector{Int}(undef, 2)
#     rowvals2 = Vector{Int}(undef, L2)
#     v2 = SparseMatrixCSC(length(gridmask), 1, colptr2, rowvals2, nzvals2)
#
#     gridmaskcache2 = falses(size(b̃b̃))
#     rowvalscache2 = Vector{Int}(undef,L2)
#     nzvalscache2 = Vector{T}(undef,L2)
#     indexcache2 = Vector{Int}(undef,L2)
#
#
#     for (i,k) in enumerate(indices)
#         for (j,d) in enumerate(dff)
#             R[j,i] = - _sparse_innerproduct!(bb, b̃b̃, b_support, b̃_support, m, k, d, gridmask, gridsize,
#                 v1, gridmaskcache1, rowvalscache1, nzvalscache1,indexcache1,
#                 v2, gridmaskcache2, rowvalscache2, nzvalscache2,indexcache2)
#
#             if d==zerocartesian
#                 R[j,i] += 1
#             end
#         end
#     end
#     R
# end
#
#
# function _sparse_innerproduct(bb::Array{T,N}, b̃b̃::Array{T,N}, b_support::CartesianIndices{N},
#         b̃_support::CartesianIndices{N}, m::NTuple{N,Int}, k::CartesianIndex{N},
#         d::CartesianIndex{N}, gridmask::BitArray{N}, gridsize::NTuple{N,Int}) where {T,N}
#     L1 = length(bb)
#     nzvals1 = Vector{T}(undef, L1)
#     colptr1 = Vector{Int}(undef, 2)
#     rowvals1 = Vector{Int}(undef, L1)
#     v1 = SparseMatrixCSC(length(gridmask), 1, colptr1, rowvals1, nzvals1)
#
#     gridmaskcache1 = falses(size(bb))
#     rowvalscache1 = Vector{Int}(undef,L1)
#     nzvalscache1 = Vector{T}(undef,L1)
#     indexcache1 = Vector{Int}(undef,L1)
#
#     L2 = length(b̃b̃)
#     nzvals2 = Vector{T}(undef, L2)
#     colptr2 = Vector{Int}(undef, 2)
#     rowvals2 = Vector{Int}(undef, L2)
#     v2 = SparseMatrixCSC(length(gridmask), 1, colptr2, rowvals2, nzvals2)
#
#     gridmaskcache2 = falses(size(b̃b̃))
#     rowvalscache2 = Vector{Int}(undef,L2)
#     nzvalscache2 = Vector{T}(undef,L2)
#     indexcache2 = Vector{Int}(undef,L2)
#
#     _sparse_innerproduct!(bb, b̃b̃, b_support, b̃_support, m, d, gridmask, gridsize,
#         v1, gridmaskcache1, rowvalscache1, nzvalscache1,indexcache1,
#         v2, gridmaskcache2, rowvalscache2, nzvalscache2,indexcache2)
# end
#
# function _sparse_innerproduct!(bb::Array{T,N}, b̃b̃::Array{T,N}, b_support::CartesianIndices{N},
#         b̃_support::CartesianIndices{N}, m::NTuple{N,Int}, k::CartesianIndex{N},
#         d::CartesianIndex{N}, gridmask::BitArray{N}, gridsize::NTuple{N,Int},
#         v1, gridmaskcache1::BitArray{N}, rowvalscache1, nzvalscache1, indexcache1,
#         v2, gridmaskcache2::BitArray{N}, rowvalscache2, nzvalscache2, indexcache2) where {T,N}
#
#     sparsevector!(v1, b_support .+ CartesianIndex((k.I .-1 ) .* m), gridsize, gridmask, bb, gridmaskcache1, rowvalscache1, nzvalscache1,indexcache1)
#     sparsevector!(v2, b̃_support .+ CartesianIndex(m .* (d.I .+ k.I .- 1)), gridsize, gridmask, b̃b̃, gridmaskcache2, rowvalscache2, nzvalscache2,indexcache2)
#
#     dot(v1,v2)
# end
#
# function sparsevector!(v, indices::CartesianIndices{N}, gridsize::NTuple{N,Int}, gridmask::BitArray{N}, bb::Array{T,N}, gridmaskcache::BitArray{N},
#         rowvalscache, nzvalscache, indexcache) where {T,N}
#     mod_indices = ModCartesianIndices(gridsize, first(indices), last(indices))
#     L = LinearIndices(size(gridmask))
#     cacheindices = CartesianIndices(size(mod_indices))
#
#     copyto!(gridmaskcache, cacheindices, gridmask,mod_indices)
#     v.colptr[1] = 1
#
#     nzval = 0
#     rowval = 0
#
#     for (i,j) in zip(cacheindices,mod_indices)
#         rowval += 1
#         if gridmaskcache[i]
#             nzval += 1
#             rowvalscache[nzval] = L[j]
#             nzvalscache[nzval] = bb[i]
#         end
#     end
#     v.colptr[2] = nzval + 1
#
#     for i in nzval+1:length(nzvalscache)
#         rowvalscache[i] = length(L)+1
#     end
#     sortperm!(indexcache, rowvalscache)
#     # sortperm!(view(indexcache,1:nzval), view(rowvalscache,1:nzval))
#     for i in 1:nzval
#         v.rowval[i] = rowvalscache[indexcache[i]]
#         v.nzval[i] = nzvalscache[indexcache[i]]
#     end
#     v
# end
#
# splitted(M::ModCartesianIndices{N}) where N =
#     first(M) > last(M)
#
# function Base.copyto!(dest::AbstractArray{T,N}, destindices::CartesianIndices, src::AbstractArray{T,N}, srcindices::ModCartesianIndices) where {T,N}
#     if !splitted(srcindices)
#         return copyto!(dest, destindices, src, first(srcindices):last(srcindices))
#     end
#
#     p1 = first(srcindices)
#     p2 = CartesianIndex(srcindices.size)
#     p3 = last(srcindices)
#
#     srcix = p1:p2
#
#     q1 = first(destindices)
#     q2 = q1 + CartesianIndex(size(srcix) .- 1)
#     q3 = last(destindices)
#
#     destix = q1:q2
#     copyto!(dest, destix, src, srcix)
#
#     p2 = p2 + CartesianIndex{N}(1)
#     q2 = q2 + CartesianIndex{N}(1)
#
#     for i in CartesianIndices(ntuple(k->2,Val(N)))
#         if i == CartesianIndex{N}(1)
#             nothing
#         else
#             srcix = CartesianIndex(ntuple(k->i[k] == 1 ? p1[k] : p2[k] ,Val(N))):CartesianIndex(ntuple(k->i[k] == 1 ? p2[k] : p3[k],Val(N)))
#             destix = CartesianIndex(ntuple(k->i[k] == 1 ? q1[k] : q2[k] ,Val(N))):CartesianIndex(ntuple(k->i[k] == 1 ? q2[k] : q3[k],Val(N)))
#             copyto!(dest, destix, src, srcix)
#         end
#
#     end
#
#     dest
# end
#
# function sparsemixedgramcomplement(indices::AbstractArray{CartesianIndex{N}},
#         b̃::NTuple{N,CompactInfiniteVector}, b::NTuple{N,CompactInfiniteVector}, m::NTuple{N,Int}, grid::AbstractGrid;atol=1e-14) where N
#     T = promote_type(map(eltype,b)..., map(eltype, b̃)...)
#     b_support = support(b)
#     b̃_support = support(b̃)
#     dff = difference_indices(b,b̃,m)
#
#     b_supportlength = length(dff)
#
#     R = sparsemixedgramcomplement_nzband(indices, T, b̃, b̃_support, b, b_support, m, dff, grid)
#
#     nonzeromaskR = abs.(R) .> atol
#     nnz = count(nonzeromaskR)
#     basissize = div.(size(supergrid(grid)),m)
#     n = length(indices)
#     m = prod(basissize)
#     nzvals = Vector{T}(undef, nnz)
#     colptr = Vector{Int}(undef, n+1)
#     rowvals = Vector{Int}(undef, nnz)
#
#     rowvalscol = Vector{Int}(undef, length(dff))
#     nzvalscol = Vector{T}(undef, length(dff))
#     colix = Vector{Int}(undef, length(dff))
#
#     L = LinearIndices(basissize)
#     # periodic = ntuple(k->true,Val(N))
#     valn = 1
#     colptr[1] = 1
#     @inbounds for (i,k) in enumerate(indices)
#         nnzcol = 0
#         kdff = k .+ dff
#         ls = ModCartesianIndices(basissize,first(kdff),last(kdff))
#
#         for (j,l) in enumerate(ls)
#             if nonzeromaskR[j,i]
#                 nnzcol += 1
#                 nzvalscol[nnzcol] = R[j,i]
#                 rowvalscol[nnzcol] = L[l]
#             end
#         end
#         for j in nnzcol+1:length(dff)
#             rowvalscol[j] = length(L)+1
#         end
#         for j in 1:b_supportlength
#             colix[j] = j
#         end
#         sort!(colix ,1,b_supportlength, InsertionSort,Base.Perm(Base.Order.ForwardOrdering(),rowvalscol))
#         # sortperm!(colix,rowvalscol)
#         for j in 1:nnzcol
#             nzvals[valn] = nzvalscol[colix[j]]
#             rowvals[valn] = rowvalscol[colix[j]]
#             valn += 1
#         end
#
#         colptr[i+1] = colptr[i] + nnzcol
#     end
#     SparseMatrixCSC(m,n,colptr,rowvals,nzvals)
# end
#
# @inline function _innerproduct(bb::Array{T,N}, b̃b̃::Array{T,N}, b_support::CartesianIndices{N},
#     b̃_support::CartesianIndices{N},
#     m::NTuple{N,Int}, k::CartesianIndex{N}, coef_dff::CartesianIndex{N}, gridmask::BitArray{N}, gridsize::NTuple{N,Int}) where {N,T}
#
#     grid_dff = CartesianIndex(m .* coef_dff.I)
#     overlapping_support = overlapping(first(b̃_support) + grid_dff, last(b̃_support) + grid_dff,first(b_support),last(b_support))
#
#     ktranslate = CartesianIndex((k.I .-1 ) .* m)
#
#     overlapping_k_support = overlapping_support .+ ktranslate
#     # prepare to put overlapping_k_support as index in maskedgrid
#     I = GridArrays.ModCartesianIndicesBase.ModCartesianIndices(gridsize, first(overlapping_k_support), last(overlapping_k_support))
#
#     r = T(0)
#     @inbounds for (i,imod) in zip(overlapping_k_support,I)
#         if gridmask[imod]
#             k = i - ktranslate
#             if k∈b_support
#                 b1 = bb[(k -first(b_support)+CartesianIndex{N}(1))]
#                 k = i-ktranslate-grid_dff
#                 if k∈b̃_support
#                     b2 = b̃b̃[(k -first(b̃_support)+CartesianIndex{N}(1))]
#                     r += b1*b2
#                 end
#             end
#         end
#     end
#     r
# end
#
# function overlapping(a::CartesianIndices, b::CartesianIndices)
#     l1,r1 = first(a), last(a)
#     l2,r2 = first(b), last(b)
#     overlapping(l1,r1,l2,r2)
# end
#
# function overlapping(l1,r1,l2,r2)
#     if l1 <= l2 <= r1
#         return l2:r1
#     elseif l2 <= l1 <= r2
#         return l1:r2
#     else
#         return r1:l1
#     end
#
# end
