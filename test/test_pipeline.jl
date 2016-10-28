# Tests for pipeline.jl


test_i1 = testfile("1CLL_H.pdb")
test_d1 = testfile("1CLL.dssp")
test_i2 = testfile("1CTR_H.pdb")
test_d2 = testfile("1CTR.dssp")
test_extra_pdb_1 = testfile("1CFF_1.pdb")
test_extra_pdb_2 = testfile("1CFF_2.pdb")
test_pocket_points = testfile("1CTR_pocket_points.pdb")
test_tmscore_path = "TMscore"


@testset "Pipeline" begin
    # Run the whole pipeline as an integration test
    # Temp directory which is removed at the end
    temp_dir = "$(tempdir())/$(defaults["out_dir"])"
    n_strucs = 4
    n_mods = 1

    # Two structures
    runpipeline(
        i1=test_i1,
        d1=test_d1,
        i2=test_i2,
        d2=test_d2,
        out_dir=temp_dir,
        n_strucs=n_strucs,
        extra_pdbs=[test_extra_pdb_1, test_extra_pdb_2],
    )
    # Check that output files are produced
    @test length(readdir("$temp_dir/pdbs")) == n_strucs
    rm(temp_dir, recursive=true)

    # One structure
    runpipeline(
        i1=test_i1,
        d1=test_d1,
        out_dir=temp_dir,
        n_strucs=n_strucs,
        extra_pdbs=[test_extra_pdb_1, test_extra_pdb_2],
    )
    @test length(readdir("$temp_dir/pdbs")) == n_strucs
    rm(temp_dir, recursive=true)

    # Two structures with perturbation
    runpipeline(
        i1=test_i1,
        d1=test_d1,
        i2=test_i2,
        d2=test_d2,
        out_dir=temp_dir,
        n_strucs=n_strucs,
        extra_pdbs=[test_extra_pdb_1, test_extra_pdb_2],
        mod_path=test_pocket_points,
        n_mods=n_mods,
    )
    @test length(readdir("$temp_dir/pdbs_mod_1")) == n_strucs
    rm(temp_dir, recursive=true)

    # One structure with perturbation
    runpipeline(
        i1=test_i2,
        d1=test_d2,
        out_dir=temp_dir,
        n_strucs=n_strucs,
        extra_pdbs=[test_extra_pdb_1, test_extra_pdb_2],
        mod_path=test_pocket_points,
        n_mods=n_mods,
    )
    @test length(readdir("$temp_dir/pdbs_mod_1")) == n_strucs
    rm(temp_dir, recursive=true)

    @test_throws AssertionError runpipeline(i2=test_i2, d2=test_d2)

    # Auto-parameterisation pipeline
    # Answer non-deterministic so cannot check it directly
    parampipeline(
        i1=test_i1,
        d1=test_d1,
        i2=test_i2,
        d2=test_d2,
        out_dir=temp_dir,
        n_strucs=n_strucs,
        tmscore_path=test_tmscore_path,
    )
    @test length(readdir("$temp_dir/tw_1_0/pdbs")) == n_strucs
    rm(temp_dir, recursive=true)

    @test_throws AssertionError parampipeline(i1=test_i1, d1=test_d1)
    @test_throws AssertionError parampipeline(i1=test_i1, d1=test_d1,
        i2=test_i2, d2=test_d2, tmscore_path="invalidpath")
end
