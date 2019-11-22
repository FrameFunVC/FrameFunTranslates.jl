module FrameFunTranslates

using Reexport

module CompactInfiniteVectors
    import CompactTranslatesDict.CompactInfiniteVectors: compactinfinitevector
    using FrameFun, GridArrays
    compactinfinitevector(dict::ExtensionFrame, grid::AbstractSubGrid) = compactinfinitevector(basis(dict), supergrid(grid))
end

include("SPQR_Solvers.jl")
@reexport using .SPQR_Solvers

include("SparseArrayOperators.jl")
@reexport using .SparseArrayOperators

include("TranslatesPlatforms.jl/TranslatesPlatforms.jl")
@reexport using .TranslatesPlatforms

include("CompactAZ.jl/CompactAZ.jl")
@reexport using .CompactAZ

end
