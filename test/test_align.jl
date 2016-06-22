# Tests for align.jl


@testset "Align" begin
    coords_one = [
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 0.0;
    ]
    coords_two = [
        0.0 -1.0 0.0;
        1.0 0.0 0.0;
        1.0 1.0 1.0;
    ]
    trans_one_test, trans_two_test, rotation_test = kabschalignment(coords_one, coords_two)
    trans_one_real = [1/3, 1/3, 0]
    trans_two_real = [-1/3, 1/3, 1]
    rotation_real = [
        0.0 -1.0 0.0;
        1.0 0.0 0.0;
        0.0 0.0 1.0;
    ]
    @test isapprox(trans_one_test, trans_one_real)
    @test isapprox(trans_two_test, trans_two_real)
    @test isapprox(rotation_test, rotation_real)


    coords_one = [
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 0.0;
    ]
    coords_two = [
        0.0 -1.0 0.0;
        1.0 0.0 0.0;
        1.0 1.0 1.0;
    ]
    devs_test = displacements(coords_one, coords_two)
    devs_real = [sqrt(3), sqrt(3), 1.0]
    @test isapprox(devs_test, devs_real)


    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_thr = Atom("CA", "ALA", 'A', 21, [0.0, 0.0, 0.0], "C")
    atoms = [atom_one, atom_two, atom_thr]
    coords_one = [
        0.0 1.0 2.0;
        0.0 1.0 0.0;
        0.0 0.0 0.0;
    ]
    coords_two = [
        0.0 -1.0 0.0;
        0.0 0.0 3.0;
        3.0 3.0 3.0;
    ]
    align!(coords_one, coords_two, atoms)
    coords_one_real = [
        0.0 0.5 3.0;
        -1.0 1.5 3.0;
        0.0 2.5 3.0;
    ]
    coords_one_real = [
        0.0 -1.0 0.0;
        0.5 1.5 2.5;
        3.0 3.0 3.0;
    ]
    @test coords_one == coords_one_real


    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_thr = Atom("CA", "ALA", 'A', 21, [0.0, 0.0, 0.0], "C")
    atoms = [atom_one, atom_two, atom_thr]
    coords_one = [
        0.0 1.0 2.0;
        0.0 1.0 0.0;
        0.0 0.0 0.0;
    ]
    coords_two = [
        0.0 -1.0 0.0;
        0.0 0.0 3.0;
        3.0 3.0 3.0;
    ]
    alignsimple!(coords_one, coords_two, atoms)
    coords_one_real = [
        0.0 -1.0 0.0;
        0.5 1.5 2.5;
        3.0 3.0 3.0;
    ]
    @test coords_one == coords_one_real


    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("C", "ALA", 'A', 20, [1.0, 1.0, 0.0], "C")
    atom_thr = Atom("CA", "ALA", 'A', 21, [2.0, 0.0, 0.0], "C")
    atoms = [atom_one, atom_two, atom_thr]
    atoms_new = deepcopy(atoms)
    coords_ref = [
        0.0 -1.0 0.0;
        0.0 0.0 3.0;
        3.0 3.0 3.0;
    ]
    alignatoms!(atoms, coords_ref)
    coords_real = [
        0.0 -1.0 0.0;
        0.5 1.5 2.5;
        3.0 3.0 3.0;
    ]
    @test atomcoords(atoms) == coords_real
    atoms_ref = atomcoords(atoms, coords_ref)
    alignatoms!(atoms_new, atoms_ref)
    @test atomcoords(atoms_new) == coords_real


    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CA", "ALA", 'A', 30, [0.0, 0.0, 0.0], "C")
    struc_one = ModelledStructure(0.0, [0.0 5.0; 0.0 5.0; 0.0 5.0])
    struc_two = ModelledStructure(0.0, [1.0 4.0; 1.0 4.0; 1.0 4.0])
    ensemble = ModelledEnsemble([atom_one, atom_two], [struc_one, struc_two])
    average_test = centroid(ensemble)
    average_real = [0.5 4.5; 0.5 4.5; 0.5 4.5]
    @test average_test == average_real


    coords_one = [
        0.0 3.0 13.0;
        0.0 0.0 0.0;
        0.0 0.0 0.0;
    ]
    coords_two = [
        0.0 5.0 10.0;
        1.0 0.0 0.0;
        0.0 0.0 0.0;
    ]
    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CA", "ALA", 'A', 30, [0.0, 0.0, 0.0], "C")
    atom_thr = Atom("C", "ALA", 'A', 40, [0.0, 0.0, 0.0], "C")
    atoms = [atom_one, atom_two, atom_thr]
    @test isapprox(rmsd(coords_one, coords_two, atoms), sqrt(2.5))


    coords_one = [
        0.0 3.0;
        0.0 0.0;
        0.0 0.0;
    ]
    coords_two = [
        0.0 5.0;
        1.0 0.0;
        0.0 0.0;
    ]
    @test isapprox(rmsd(coords_one, coords_two), sqrt(2.5))
end
