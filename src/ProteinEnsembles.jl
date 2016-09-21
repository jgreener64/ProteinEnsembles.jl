__precompile__()

module ProteinEnsembles

if VERSION < v"0.5-"
    typealias ASCIIString String
end

include("types.jl")
include("defaults.jl")
include("bonds.jl")
include("atoms.jl")
include("io.jl")
include("interactions.jl")
include("generate.jl")
include("align.jl")
include("pca.jl")
include("perturb.jl")
include("plot.jl")
include("pipeline.jl")

end # ProteinEnsembles
