module CompactAZ

using Reexport

include("ReducedAZ.jl/ReducedAZ.jl")
@reexport using .ReducedAZ

include("SparseAZ.jl/SparseAZ.jl")
@reexport using .SparseAZ

end
