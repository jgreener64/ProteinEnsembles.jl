# Tests for pipeline.jl


test_i1 = testfile("4AKE_H.pdb")
test_d1 = testfile("4AKE.dssp")
test_i2 = testfile("1AKE_H.pdb")
test_d2 = testfile("1AKE.dssp")
test_extra_pdb = testfile("1AK2.pdb")
test_pocket_points = testfile("1AKE_pocket_points.pdb")


@testset "Pipeline" begin
    # Run the whole pipeline as an integration test
    # Temp directory which is removed at the end
    temp_dir = "$(tempdir())/$(defaults["out_dir"])"

    # Two structures
    runpipeline(
        i1=test_i1,
        d1=test_d1,
        i2=test_i2,
        d2=test_d2,
        out_dir=temp_dir,
        n_strucs=4,
        extra_pdbs=[test_extra_pdb],
    )
    # Check that output files are produced
    @test length(readdir("$temp_dir/pdbs")) == 4
    rm(temp_dir, recursive=true)

    # One structure
    runpipeline(
        i1=test_i1,
        d1=test_d1,
        out_dir=temp_dir,
        n_strucs=4,
        extra_pdbs=[test_extra_pdb],
    )
    @test length(readdir("$temp_dir/pdbs")) == 4
    rm(temp_dir, recursive=true)

    # Two structures with perturbation
    runpipeline(
        i1=test_i1,
        d1=test_d1,
        i2=test_i2,
        d2=test_d2,
        out_dir=temp_dir,
        n_strucs=4,
        extra_pdbs=[test_extra_pdb],
        mod_path=test_pocket_points,
        n_mods=2,
    )
    @test length(readdir("$temp_dir/pdbs_mod_1")) == 4
    rm(temp_dir, recursive=true)

    # One structure with perturbation
    runpipeline(
        i1=test_i2,
        d1=test_d2,
        out_dir=temp_dir,
        n_strucs=4,
        extra_pdbs=[test_extra_pdb],
        mod_path=test_pocket_points,
        n_mods=2,
    )
    @test length(readdir("$temp_dir/pdbs_mod_1")) == 4
    rm(temp_dir, recursive=true)

    @test_throws AssertionError runpipeline(i2=test_i2, d2=test_d2)
end
