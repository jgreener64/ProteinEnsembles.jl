using ProteinEnsembles
using BaseTestNext

const Test = BaseTestNext

println("These tests usually take about 5 minutes to run")

testfile(path::AbstractString...) = Pkg.dir("ProteinEnsembles", "test", "test_files", path...)

include("test_types.jl")
include("test_atoms.jl")
include("test_io.jl")
include("test_interactions.jl")
include("test_generate.jl")
include("test_align.jl")
include("test_pca.jl")
include("test_perturb.jl")
include("test_pipeline.jl")
