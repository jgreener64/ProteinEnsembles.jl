using ProteinEnsembles
using Test


testfile(path::AbstractString...) = joinpath(dirname(@__FILE__), "test_files", path...)


# Whether to run the parameterisation tests
# This is to satisfy the online auto-builds
# Set to true locally to run the tests
run_param_test = false


if run_param_test
    # Optional argument is path to TMscore executable
    if isempty(ARGS)
        test_tmscore_path = "TMscore"
    else
        test_tmscore_path = ARGS[1]
    end
    @assert tmscorepathvalid(test_tmscore_path) "Not a valid TMscore path: \"$test_tmscore_path\""
    println("TMscore executable path is taken as \"$test_tmscore_path\"")
end


# Temp directory which is removed at the end
temp_dir = "$(tempdir())/$(defaults["out_dir"])"


include("test_types.jl")
include("test_atoms.jl")
include("test_io.jl")
include("test_interactions.jl")
include("test_generate.jl")
include("test_align.jl")
include("test_pca.jl")
include("test_perturb.jl")
include("test_pipeline.jl")
