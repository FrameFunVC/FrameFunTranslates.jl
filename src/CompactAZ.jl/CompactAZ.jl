module CompactAZ

using Reexport

module CompactFrameFunExtension
    using Reexport, FrameFun, BasisFunctions, FrameFun.ExtensionFramePlatforms, FrameFun.ApproximationProblems, InfiniteVectors, SparseArrays, ....TranslatesPlatforms
    using FrameFun.FrameFunInterface: @trial, @aptoplatform
    using ....SPQR_Solvers

    macro efplatformtobasisplatforms(ex)
        efex = Meta.parse("ef_"*string(ex))
        ret = quote
            $(ex)(ss::SamplingStyle, platform::ExtensionFramePlatform, param, args...; options...) =
                $(efex)(ss, platform, param, platform.basisplatform, args...; options...)
            $(efex)(ss::SamplingStyle, platform::ExtensionFramePlatform, param, bplatform::ProductPlatform, args...; options...) =
                $(efex)(ss, platform, param, elements(bplatform), args...; options...)
            $(efex)(ss::SamplingStyle, platform::ExtensionFramePlatform, param, bplatform::BasisPlatform, args...; options...) =
                $(efex)(ss, platform, param, tuple(bplatform), args...; options...)
        end
        esc(ret)
    end

    export CompactAZStyle
    abstract type CompactAZStyle <: SolverStyle end
    export SparseAZStyle
    """
        struct SparseAZStyle <: CompactAZStyle end

    Use the sparse AZ algorithm to solve the discretized approximation problem.
    """
    struct SparseAZStyle <: CompactAZStyle end

    export ReducedAZStyle
    """
        struct ReducedAZStyle <: CompactAZStyle end

    Use the reduced AZ algorithm to solve the discretized approximation problem.
    """
    struct ReducedAZStyle <: CompactAZStyle end


    export reducedAAZAoperator
    @trial reducedAAZAoperator
    @efplatformtobasisplatforms reducedAAZAoperator
    reducedAAZAoperator(samplingstyle::SamplingStyle, ap::ApproximationProblem; L=samplingparameter(ap), options...) =
        reducedAAZAoperator(samplingstyle, ap, L; options...)
    ef_reducedAAZAoperator(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L; options...) =
        default_ef_reducedAAZAoperator(samplingstyle, platform, param, platforms, L; return_nzrows=false, options...)

    function default_ef_reducedAAZAoperator(samplingstyle::SamplingStyle, platform::Platform, param, platforms, L; verbose=false, return_nzrows, options...)
        nonzero_coefs = haskey(options, :nonzero_coefs) ? options[:nonzero_coefs] : ef_nonzero_coefficients(samplingstyle, platform, param, platforms, L; verbose=verbose, options...)
        M = haskey(options, :AAZA) ? options[:AAZA] : firstAZstepoperator(platform, param; samplingstyle=samplingstyle, L=L, verbose=verbose, options...)
        nz_tol = haskey(options, :nz_tol) ? options[:nz_tol] : default_threshold(M)

        verbose && @info "Reduced AZ: Restrict columns (coefficients) from $(size(src(M))) to $(length(nonzero_coefs))"
        s = zeros(src(M))
        d = zeros(dest(M))
        m = Array{eltype(M)}(undef, size(M,1), length(nonzero_coefs))

        for (j,i) in enumerate(nonzero_coefs)
            s[i] = 1
            apply!(M, d, s)
            s[i] = 0
            copyto!(m, size(M,1)*(j-1)+1, d, 1)
        end

        verbose && @info "Reduced AZ: Restricting rows..."
        nz_rows = findall(nonzero_rows(m,size(dest(M));nz_tol=nz_tol))[:]
        I = LinearIndices(size(dest(M)))[nz_rows]

        if length(nz_rows) < size(M,1)
            m = m[I,:]
            verbose && @info "Reduced AZ: Restrict rows (collocation points) from $(size(dest(M))) to $(length(nz_rows))"
        else
            verbose && @info "Reduced AZ: Rows are not restricted"
        end

        if return_nzrows
            ArrayOperator(m), nz_rows
        else
            ArrayOperator(m)
        end
    end

    export reducedAZ_AAZAreductionsolver
    @trial reducedAZ_AAZAreductionsolver
    @efplatformtobasisplatforms reducedAZ_AAZAreductionsolver
    function reducedAZ_AAZAreductionsolver(ss::SamplingStyle, ap::ApproximationProblem; L=samplingparameter(ap), REG=pQR_solver, options...)
        if haskey(options, :directsolver)
            @warn "Are you sure you want to use the `:directsolver` options and not `REG` in an AZ algorithm? "
        end
        reducedAZ_AAZAreductionsolver(ss, ap, L, REG; options...)
    end

    ef_reducedAZ_AAZAreductionsolver(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L, REG; options...) =
        default_ef_reducedAZ_AAZAreductionsolver(samplingstyle, platform, param, platforms, L, REG; options...)
    function default_ef_reducedAZ_AAZAreductionsolver(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L, REG; verbose=false, options...)
        nonzero_coefs = haskey(options, :nonzero_coefs) ? options[:nonzero_coefs] : ef_nonzero_coefficients(samplingstyle, platform, param, platforms, L;verbose=verbose,  options...)
        M = firstAZstepoperator(platform, param; samplingstyle=samplingstyle, L=L, options...)
        rM, nonzero_rows = ef_reducedAAZAoperator(samplingstyle, platform, param, platforms, L; verbose=verbose, AAZAoperator=M, nonzero_coefs=nonzero_coefs, return_nzrows=true, options...)

        dict_resop = IndexRestrictionOperator(src(M), src(rM), nonzero_coefs)
        if size(rM,1) < size(M,1)
            I = LinearIndices(size(dest(M)))[nonzero_rows]
            grid_resop = IndexRestrictionOperator(dest(M), GridBasis{coefficienttype(dest(M))}(grid(dest(M))[I]),nonzero_rows)
        else
            grid_resop = IdentityOperator(dest(M),dest(rM))
        end

        verbose && @info "Reduced AZ: use $(REG) as solver for first reduced AZ step"
        (dict_resop')*REG(rM; verbose=verbose, options...)*grid_resop
    end

    using GridArrays, CompactTranslatesDict.CompactInfiniteVectors
    using ....TranslatesPlatforms.BSplinePlatforms: AbstractPeriodicEquispacedTranslatesPlatform
    export compactsupport
    @trial compactsupport
    compactsupport(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
        compactsupport(samplingstyle, platform(ap), parameter(ap), samplingparameter(ap); options...)

    function compactsupport(ss::DiscreteStyle, platform::Platform, param, L; options...)
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L; options...)
        compactsupport(ss, basisplatform, param, os_grid; options..., os_grid=os_grid)
    end
    compactsupport(ss::DiscreteStyle, bplatform::ExtensionFramePlatform, param, platforms::Tuple, os_grid::AbstractGrid; options...) =
        error("Should not reach here.")
    compactsupport(ss::DiscreteStyle, bplatform::Platform, param, platforms::Tuple, os_grid::AbstractGrid; options...) =
        compactsupport(compactinfinitevectors(ss, bplatform, param, platforms, os_grid; options...))
    function compactsupport(vecs::NTuple{N,CompactInfiniteVector}) where N
        supports = map(InfiniteVectors.support, vecs)
        CartesianIndex(map(x->x[1], supports)):CartesianIndex(map(x->x[end], supports))
    end


    export compactinfinitevectors
    @trial compactinfinitevectors
    compactinfinitevectors(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
        compactinfinitevectors(samplingstyle, platform(ap), parameter(ap), samplingparameter(ap); options...)
    function compactinfinitevectors(ss::DiscreteStyle, platform::Platform, param, L; options...)
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L; options...)
        compactinfinitevectors(ss, platform, param, os_grid; options...)
    end
    compactinfinitevectors(ss::DiscreteStyle, bplatform::ProductPlatform, param, os_grid::AbstractGrid; options...) =
        compactinfinitevectors(ss, bplatform, param, elements(bplatform), os_grid; options...)
    compactinfinitevectors(ss::DiscreteStyle, bplatform::BasisPlatform, param, os_grid::AbstractGrid; options...) =
        compactinfinitevectors(ss, bplatform, param, tuple(bplatform), os_grid; options...)
    compactinfinitevectors(ss::DiscreteStyle, bplatform::Platform, param, platforms::Tuple{<:AbstractPeriodicEquispacedTranslatesPlatform}, os_grid::AbstractIntervalGrid; options...) =
        tuple(compactinfinitevector(dictionary(bplatform, param), os_grid; options...))
    compactinfinitevectors(ss::DiscreteStyle, bplatform::Platform, param, platforms::Tuple{Vararg{<:AbstractPeriodicEquispacedTranslatesPlatform}}, os_grid::ProductGrid; options...) =
        map((x,y)->compactinfinitevector(x,y; options...), map(dictionary, platforms, param), elements(os_grid))

    export nonzero_coefficients
    @trial nonzero_coefficients
    @efplatformtobasisplatforms nonzero_coefficients
    nonzero_coefficients(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
        nonzero_coefficients(samplingstyle, platform(ap), parameter(ap), samplingparameter(ap); options...)

    function ef_nonzero_coefficients(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L; options...)
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L; options...)
        q = div.(L, size(dictionary(platform,param)))
        ef_nonzero_coefficients(samplingstyle, platform, param, platforms, q, os_grid; options...)
    end
    ef_nonzero_coefficients(ss::DiscreteStyle, platform::Platform, param, platforms::Tuple, q, os_grid::AbstractGrid; dict=dictionary(platform,param), options...) =
        default_ef_nonzero_coefficients(ss, dict)
    function default_ef_nonzero_coefficients(ss::SamplingStyle, dict::Dictionary; options...)
        @debug "Not enough information to detect sparse structure"
        error()
        eachindex(dict)
    end

    ef_nonzero_coefficients(ss::DiscreteStyle, platform::Platform, param, platforms::Tuple{Vararg{<:AbstractPeriodicEquispacedTranslatesPlatform}}, q, os_grid::AbstractGrid; options...) =
        _nonzero_coefficients(compactsupport(ss, platform.basisplatform, param, platforms, supergrid(os_grid); os_grid=supergrid(os_grid), options...), q, mask(os_grid))

    using GridArrays.ModCartesianIndicesBase: ModCartesianIndices
    function _nonzero_coefficients(b_support::CartesianIndices{N}, m::Union{Int,NTuple{N,Int}}, gridmask::BitArray{N}) where N
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

    include("nonzero_rows.jl")

    export nonzero_pointsindices
    @trial nonzero_pointsindices
    @efplatformtobasisplatforms nonzero_pointsindices
    nonzero_pointsindices(samplingstyle::SamplingStyle, ap::ApproximationProblem, ix, relative::Bool; options...) =
        nonzero_pointsindices(samplingstyle, platform(ap), parameter(ap), samplingparameter(ap), ix, relative; options...)

    function ef_nonzero_pointsindices(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L, ix, relative::Bool; options...)
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L; options...)
        q = div.(L, size(dictionary(platform, param)))
        ef_nonzero_pointsindices(samplingstyle, platform, param, platforms, q, os_grid, ix, relative; options...)
    end
    ef_nonzero_pointsindices(ss::DiscreteStyle, platform::Platform, param, platforms::Tuple{Vararg{<:AbstractPeriodicEquispacedTranslatesPlatform}}, q, os_grid::AbstractGrid, ix, relative::Bool; options...) =
        _nonzero_pointsindices(compactsupport(ss, platform.basisplatform, param, platforms, supergrid(os_grid); os_grid=supergrid(os_grid), options...), q, mask(os_grid), ix, relative)

    function _nonzero_pointsindices(b_support::CartesianIndices{N}, m::Union{Int,NTuple{N,Int}}, gridmask::BitArray{N}, ix, relative::Bool) where N
        gridsize = size(gridmask)
        A = falses(gridsize)
        for k in ix
            l = ((k isa Int ? k : k.I) .- 1).*m
            bk_support = b_support .+ CartesianIndex(l)
            for i in ModCartesianIndices(gridsize, first(bk_support), last(bk_support))
                A[i] = gridmask[i]
            end
        end
        if relative
            findall(A[gridmask])
        else
            findall(A)
        end
    end
    export SE_K
    const SE_K = nonzero_coefficients
    const AAZA_nonzero_column_indexset = nonzero_coefficients

    export AAZA_nonzero_row_indexset
    @trial AAZA_nonzero_row_indexset
    @efplatformtobasisplatforms AAZA_nonzero_row_indexset
    AAZA_nonzero_row_indexset(ss::SamplingStyle, ap::ApproximationProblem; L=samplingparameter(ap), options...) =
        AAZA_nonzero_row_indexset(ss, ap, L; options...)

    ef_AAZA_nonzero_row_indexset(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L; options...) =
        default_ef_AAZA_nonzero_row_indexset(samplingstyle, platform, param, platforms, L; options...)
    function default_ef_AAZA_nonzero_row_indexset(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L; verbose=false, options...)
        WE_K = haskey(options, :nonzero_coefs) ? options[:nonzero_coefs] : ef_nonzero_coefficients(samplingstyle, platform, param, platforms, L;verbose=verbose,  options...)
        verbose && @info "ReducedAZStyle: WE_K has $(length(WE_K)) elements"
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L; verbose=verbose, options...)
        verbose && @info "ReducedAZStyle: grid has $(length(os_grid)) elements"
        # Oversampling
        gridmask = mask(os_grid)

        # WE_I1 is the union of two sets. The first one:
        WE_I1a = WE_K
        # The second one checks the overlapping supports of the dual
        frame2 = azdual_dict(platform, param; verbose=verbose, samplingstyle=samplingstyle, L=L, options...)
        dict2 = basis(frame2)
        cvecs_dual =  dict2 isa Dictionary1d ?
            tuple(compactinfinitevector(dict2, supergrid(os_grid); verbose=verbose, options...)) :
            map((x,y)->compactinfinitevector(x, y; verbose=verbose, options...), elements(dict2), elements(supergrid(os_grid)))
        # and primal bases
        frame1 = dictionary(platform, param)
        dict1 = basis(frame1)
        cvecs =  dict1 isa Dictionary1d ?
            tuple(compactinfinitevector(dict1, supergrid(os_grid); verbose=verbose, options...)) :
            map((x,y)->compactinfinitevector(x, y; verbose=verbose, options...), elements(dict1), elements(supergrid(os_grid)))
        cvecs_supp = compactsupport(cvecs)
        # for l in WE_K
        # The primal basis has nonzero points for l in WE_K
        N = size(dict1)
        L = size(gridmask)
        q = div.(L,N)

        nonzero_pointsindices = _nonzero_pointsindices(cvecs_supp, q, gridmask, WE_K, true)
        WE_I1b = overlappingindices(cvecs_dual, os_grid, nonzero_pointsindices, N, L)
        # The union is
        WE_I1 = union(WE_I1a,WE_I1b) # Not longer sorted
        verbose && @info "ReducedAZStyle: WE_I1 has $(length(WE_I1)) elements"

        # Get the nonzero points of this new index set
        WE_M = _nonzero_pointsindices(cvecs_supp, q, gridmask, WE_I1, true)
        verbose && @info "ReducedAZStyle: WE_M has $(length(WE_M)) elements"
        WE_M
    end

    export sparseAZ_AAZAreductionsolver
    @trial sparseAZ_AAZAreductionsolver
    @efplatformtobasisplatforms sparseAZ_AAZAreductionsolver
    function sparseAZ_AAZAreductionsolver(ss::SamplingStyle, ap::ApproximationProblem; L=samplingparameter(ap), REG=SPQR_solver, options...)
        if haskey(options, :directsolver)
            @warn "Are you sure you want to use the `:directsolver` options and not `REG` in an AZ algorithm? "
        end
        sparseAZ_AAZAreductionsolver(ss, ap, L, REG; options...)
    end

    ef_sparseAZ_AAZAreductionsolver(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L, REG; options...) =
        default_ef_sparseAZ_AAZAreductionsolver(samplingstyle, platform, param, platforms, L, REG; options...)

    function default_ef_sparseAZ_AAZAreductionsolver(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L, REG;
            verbose=false, weightedAZ=false,  options...)
        nonzero_coefs = haskey(options, :nonzero_coefs) ? options[:nonzero_coefs] : ef_nonzero_coefficients(samplingstyle, platform, param, platforms, L;verbose=verbose,  options...)
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L; verbose=verbose, options...)
        rM = haskey(options, :sparse_reducedAAZAoperator) ? options[:sparse_reducedAAZAoperator] : sparse_reducedAAZAoperator(samplingstyle, platform, param, L; verbose=verbose, os_grid=os_grid, nonzero_coefs=nonzero_coefs, options...)
        verbose && @info "Sparse AZ: use $(REG) as solver for first sparse AZ step"



        if weightedAZ
            error("not implemented")
        end
        IndexExtensionOperator(dictionary(platform, param),nonzero_coefs)*REG(rM; verbose=verbose, options...)
    end

    export sparse_reducedAAZAoperator
    @trial sparse_reducedAAZAoperator
    @efplatformtobasisplatforms sparse_reducedAAZAoperator
    function sparse_reducedAAZAoperator(ss::SamplingStyle, ap::ApproximationProblem; L=samplingparameter(ap), options...)
        if haskey(options, :directsolver)
            @warn "Are you sure you want to use the `:directsolver` options and not `REG` in an AZ algorithm? "
        end
        sparse_reducedAAZAoperator(ss, ap, L; options...)
    end

    ef_sparse_reducedAAZAoperator(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L; options...) =
        default_ef_sparse_reducedAAZAoperator(samplingstyle, platform, param, platforms, L; options...)
    function default_ef_sparse_reducedAAZAoperator(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L; verbose=false, nz_tol=0, options...)
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L;verbose=verbose, options...)
        nonzero_coefs = haskey(options, :nonzero_coefs) ? options[:nonzero_coefs] : ef_nonzero_coefficients(samplingstyle, platform, param, platforms, L; verbose=verbose, nz_tol=nz_tol, options..., os_grid=os_grid)
        cvecs = haskey(options, :cvecs) ? options[:cvecs] : compactinfinitevectors(samplingstyle, platform.basisplatform, param, L; verbose=verbose, nz_tol=nz_tol, options..., os_grid=supergrid(os_grid))
        csupps = compactsupport(cvecs)
        verbose && @info "SparseAZStyle: create sparse A-AZ^*A"
        frame2 = azdual_dict(platform, param; options..., verbose=verbose, nz_tol=nz_tol, samplingstyle=samplingstyle, L=L)
        dict2 = basis(frame2)
        cvecs_dual =  dict2 isa Dictionary1d ?
            tuple(compactinfinitevector(dict2, supergrid(os_grid); nz_tol=nz_tol, verbose=verbose, options...)) :
            map((x,y)->compactinfinitevector(x, y; nz_tol=nz_tol, verbose=verbose, options...), elements(dict2), elements(supergrid(os_grid)))

        N = size(dict2)
        q = div.(L, N)

        ix1 = nonzero_coefs
        A = _sparseRAE(cvecs, os_grid, nonzero_coefs, N; nz_tol=nz_tol, verbose=verbose, options...)

        ix2 = overlappingindices(cvecs_dual, os_grid, findall(nonzero_rows(A)), N, L)
        ix3 = unique(sort(vcat(ix1,ix2)))

        Z = _sparseRAE(cvecs_dual, os_grid, ix3, N; nz_tol=nz_tol, verbose=verbose, options...)
        ImZA = sparseidentity(ix1,ix3)-Z'A
        droptol!(ImZA, nz_tol)
        RAE = _sparseRAE(cvecs, os_grid, ix3, N; nz_tol=nz_tol, verbose=verbose, options...)
        verbose && @info "SparseAZStyle: removing everything smaller than $nz_tol"
        M = droptol!(RAE*ImZA, nz_tol)
        verbose && @info "SparseAZStyle: A-AZ^*A has size $(size(M)) and $(nnz(M)) nonzero elements ($(100nnz(M)/prod(size(M)))% fill)"
        src = dictionary(platform, param)[nonzero_coefs]
        dest = GridBasis{coefficienttype(src)}(os_grid)
        ArrayOperator(M, src, dest)
    end

    function overlappingindices(vecs::NTuple{N,CompactInfiniteVector}, grid::AbstractGrid, ix, param, L) where N
        @assert prod(L) == prod(size(supergrid(grid)))
        boundarygrid_subindices = subindices(grid)[ix]
        if first(boundarygrid_subindices) isa Int
            boundarygrid_subindices = map(CartesianIndex, boundarygrid_subindices)
        end
        dict2mask = falses(param)
        dict2support = compactsupport(vecs)
        m = div.(L, param)
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

    @trial sparseRAE
    @efplatformtobasisplatforms sparseRAE
    sparseRAE(platform::Platform, param, ix; options...) = sparseRAE(approximationproblem(platform, param), ix; options...)
    sparseRAE(samplingstyle::SamplingStyle, ap::ApproximationProblem, ix; L=samplingparameter(ap), options...) =
        sparseRAE(samplingstyle, ap, L, ix; options...)
    function ef_sparseRAE(samplingstyle::DiscreteStyle, platform::Platform, param, platforms::Tuple, L, ix; options...)
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L; options...)
        cvecs = compactinfinitevectors(samplingstyle, platform.basisplatform, param, platforms, supergrid(os_grid); options...)
        _sparseRAE(cvecs, os_grid, ix, param; options...)
    end

    _sparseRAE(b::NTuple{1,CompactInfiniteVector}, grid::AbstractGrid, indices::AbstractVector{Int}, param::Tuple{Int}; options...) =
        _sparseRAE(b, grid, map(x->CartesianIndex{1}(x), indices), param; options...)

    _sparseRAE(b::NTuple{1,CompactInfiniteVector}, grid::AbstractGrid, indices::AbstractVector{Int}, param::Int; options...) =
        _sparseRAE(b, grid, map(x->CartesianIndex{1}(x), indices), tuple(param); options...)


    function _sparseRAE(b::NTuple{N,CompactInfiniteVector},
        grid::AbstractGrid, indices::AbstractVector{CartesianIndex{N}},
        param::NTuple{N,Int}; nz_tol=0, options...) where N
        gridsize = size(supergrid(grid))

        B = Array(OuterProductArray(map(subvector, b)...))[:]
        B[abs.(B).<nz_tol].=0
        b_support = compactsupport(b)
        b_supportlength = length(b_support)

        m = div.(gridsize,param)
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


    function ef_true_nonzero_reducedAZ_AAZAreductionsolver(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L, REG;
            verbose=false, weightedAZ=false, options...)
        if haskey(options, :directsolver)
            @warn "Are you sure you want to use the `:directsolver` options and not `REG` in an AZ algorithm? "
        end
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L; verbose=verbose,options...)
        nonzero_coefs = haskey(options, :nonzero_coefs) ? options[:nonzero_coefs] : ef_nonzero_coefficients(samplingstyle, platform, param, platforms, L; verbose=verbose, options..., os_grid=os_grid)
        rel_nonzero_points = ef_AAZA_nonzero_row_indexset(samplingstyle, platform, param, platforms, L; verbose=verbose, options..., os_grid=os_grid, nonzero_coefs=nonzero_coefs)

        rM = ef_true_nonzero_reducedAAZAoperator(samplingstyle, platform, param, platforms, L; verbose=verbose, os_grid=os_grid, nonzero_coefs=nonzero_coefs, rel_nonzero_points=rel_nonzero_points, options...)
        linop = BasisFunctions.LinearizationOperator(GridBasis(os_grid))
        dest_lin = dest(linop)

        dict = dictionary(platform, param)
        E = IndexExtensionOperator(dict,nonzero_coefs)

        if weightedAZ
            verbose && @info "Weighted AZ"
            AZ_Cweight = haskey(options,:AZ_Cweight) ? options[:AZ_Cweight] : error("No options `AZ_Cweight`")
            @assert size(src(AZ_Cweight)) == size(dict)
            W = E'*AZ_Cweight*E

            AZ_Cweight*E*
                REG(rM*W; verbose=verbose, options...)*
                    IndexRestrictionOperator(dest_lin, rel_nonzero_points)*linop
        else
            IndexExtensionOperator(dict,nonzero_coefs)*
                REG(rM; verbose=verbose, options...)*
                    IndexRestrictionOperator(dest_lin, rel_nonzero_points)*linop
        end
    end

    export true_nonzero_reducedAAZAoperator
    @trial true_nonzero_reducedAAZAoperator
    @efplatformtobasisplatforms true_nonzero_reducedAAZAoperator
    true_nonzero_reducedAAZAoperator(samplingstyle::SamplingStyle, ap::ApproximationProblem; L=samplingparameter(ap), options...) =
        true_nonzero_reducedAAZAoperator(samplingstyle, ap, L; options...)
    function ef_true_nonzero_reducedAAZAoperator(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple, L; verbose=false, options...)
        verbose && @info "Reduced AZ: Create R(A-AZA)E operator"
        os_grid = haskey(options, :os_grid) ? options[:os_grid] : sampling_grid(samplingstyle, platform, param, L;verbose=verbose, options...)
        M = firstAZstepoperator(platform, param; samplingstyle=samplingstyle, L=L, verbose=verbose, options...,)

        nonzero_coefs = haskey(options, :nonzero_coefs) ? options[:nonzero_coefs] : ef_nonzero_coefficients(samplingstyle, platform, param, platforms, L; verbose=verbose,  options..., os_grid=os_grid)
        verbose && @info "Reduced AZ: Restrict columns (coefficients) from $(length(src(M))) to $(length(nonzero_coefs))"
        WE_M = haskey(options, :rel_nonzero_points) ? options[:rel_nonzero_points] : ef_AAZA_nonzero_row_indexset(samplingstyle, platform, param, platforms, L; verbose=verbose, options..., os_grid=os_grid, nonzero_coefs=nonzero_coefs)
        verbose && @info "Reduced AZ: Restrict rows (collocation points) from $(length(os_grid)) to $(length(WE_M))"

        linop = BasisFunctions.LinearizationOperator(dest(M))
        dest_lin = dest(linop)

        grid_res = IndexRestrictionOperator(dest_lin, WE_M)*linop
        coef_res = IndexRestrictionOperator(src(M), nonzero_coefs)

        grid_res*M*coef_res'
    end


    function ef_reducedAAZAoperator(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple{Vararg{<:CDBSplinePlatform}}, L; true_nonzero=true, options...)
        if true_nonzero
            ef_true_nonzero_reducedAAZAoperator(samplingstyle, platform, param, platforms, L; options...)
        else
            default_ef_reducedAAZAoperator(samplingstyle, platform, param, platforms, L; return_nzrows=false, options...)
        end
    end

    function ef_reducedAZ_AAZAreductionsolver(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple{Vararg{<:CDBSplinePlatform}}, L, REG; true_nonzero=true, options...)
        if true_nonzero
            ef_true_nonzero_reducedAZ_AAZAreductionsolver(samplingstyle, platform, param, platforms, L, REG; options...)
        else
            default_ef_reducedAZ_AAZAreductionsolver(samplingstyle, platform, param, platforms, L, REG; options...)
        end
    end
end
@reexport using .CompactFrameFunExtension

include("ReducedAZ.jl/ReducedAZ.jl")
@reexport using .ReducedAZ

include("SparseAZ.jl/SparseAZ.jl")
@reexport using .SparseAZ

end
