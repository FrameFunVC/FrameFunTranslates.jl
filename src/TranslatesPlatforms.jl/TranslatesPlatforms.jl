module TranslatesPlatforms

using Reexport

include("BSplinePlatforms.jl")
@reexport using .BSplinePlatforms

include("NdBSplinePlatforms.jl")
@reexport using .NdBSplinePlatforms

end
