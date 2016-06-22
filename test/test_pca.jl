# Tests for pca.jl


test_ens_prefix = testfile("AKE_ensemble", "AKE.pdb")


@testset "PCA" begin
    # Values tested against Bio3D values in R
    ensemble = readensemble(test_ens_prefix, 10) # These structures are already aligned
    pca = PCA(ensemble)
    @test isapprox(pca.evals[1] / pca.evals[2], 5.26960, rtol=0.001)
    @test isapprox(abs(pca.evecs[1, 1]), 0.01500, rtol=0.001)
    @test isapprox(abs(pca.evecs[642, 1]), 0.00201, rtol=0.001)
    @test isapprox(abs(pca.evecs[1, 2]), 0.00381, rtol=0.001)


    ensemble = readensemble(test_ens_prefix, 10) # These structures are already aligned
    pca = PCA(ensemble)
    pcs = projectensemble(ensemble, pca)
    @test isapprox(abs(pcs[1, 1]), 20.32696, rtol=0.01)
    @test isapprox(abs(pcs[9, 1]), 3.73002, rtol=0.01)
    @test isapprox(abs(pcs[1, 10]), 41.06011, rtol=0.01)


    ensemble = readensemble(test_ens_prefix, 10) # These structures are already aligned
    pca = PCA(ensemble)
    pcs = projectstructure(ensemble.atoms, pca)
    @test isapprox(abs(pcs[1]), 20.32696, rtol=0.01)
    @test isapprox(abs(pcs[9]), 3.73002, rtol=0.01)


    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CA", "ALA", 'A', 21, [0.0, 0.0, 0.0], "C")
    atom_thr = Atom("CA", "ALA", 'A', 22, [0.0, 0.0, 0.0], "C")
    atoms = [atom_one, atom_two, atom_thr]
    coords_one = [
        0.0 -2.0 0.0;
        0.0 0.0 0.0;
        0.0 0.0 0.0;
    ]
    coords_two = [
        1.0 0.0 3.0;
        3.0 0.0 sqrt(3);
        0.0 0.0 0.0;
    ]
    coords_thr = [
        2.0 2.0 3.0;
        0.0 0.0 -sqrt(3);
        0.0 0.0 0.0;
    ]
    struc_one = ModelledStructure(0.0, coords_one)
    struc_two = ModelledStructure(0.0, coords_two)
    struc_thr = ModelledStructure(0.0, coords_thr)
    strucs = [struc_one, struc_two, struc_thr]
    ensemble = ModelledEnsemble(atoms, strucs)
    flucs_test = fluctuations(ensemble)
    flucs_real = [8/3, 8/3, 4.0]
    @test isapprox(flucs_test, flucs_real)
end
