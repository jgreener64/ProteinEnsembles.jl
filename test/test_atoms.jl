# Tests for atoms.jl


test_pdb_filepath = testfile("1CTR_H.pdb")
test_dssp_filepath = testfile("1CTR.dssp")


@testset "Atoms" begin
    @test inferelement("CA") == "C"
    @test inferelement("NC") == "C"
    @test inferelement("NH") == "N"
    @test inferelement("X") == "-"


    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CA", "XXX", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_three = Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    @test areatomssame(atom_one, atom_two)
    @test !areatomssame(atom_one, atom_three)


    atoms_one = Atom[
        Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C"),
        Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    ]
    atoms_two = Atom[Atom("CA", "SER", 'A', 20, [0.0, 0.0, 0.0], "C")]
    common_one_test, common_two_test = findcommonatoms(atoms_one, atoms_two)
    common_one_real = Atom[Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")]
    common_two_real = Atom[Atom("CA", "SER", 'A', 20, [0.0, 0.0, 0.0], "C")]
    @test length(common_one_test) == 1
    @test length(common_two_test) == 1
    @test common_one_test[1].atom_name == common_one_real[1].atom_name
    @test common_one_test[1].res_name == common_one_real[1].res_name
    @test common_two_test[1].atom_name == common_two_real[1].atom_name
    @test common_two_test[1].res_name == common_two_real[1].res_name


    atoms_one = [
        Atom("CA", "ALA", 'A', 22, [0.0, 0.0, 0.0], "C"),
        Atom("CA", "ALA", 'A', 19, [0.0, 0.0, 0.0], "C"),
        Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C"),
    ]
    atoms_two = [
        Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C"),
        Atom("CA", "ALA", 'A', 21, [0.0, 0.0, 0.0], "C"),
        Atom("CA", "ALA", 'A', 22, [0.0, 0.0, 0.0], "C"),
    ]
    atom_map = atommap(atoms_one, atoms_two)
    @test atom_map == [3, 0, 1]


    @test atomid(Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")) == "CA/ALA/20/A"


    atoms = readpdb(test_pdb_filepath)
    dssps = readdssp(test_dssp_filepath, atoms)
    @test length(dssps) == 142
    @test dssps["75A"] == '-'


    atom_one = Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_thr = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atoms = [atom_one, atom_two, atom_thr]
    ca_inds = calphaindices(atoms)
    @test ca_inds == [2]


    atom_one = Atom("C", "ALA", 'A', 20, [1.0, 4.0, 5.0], "C")
    atom_two = Atom("CA", "ALA", 'A', 20, [8.0, 4.0, 2.0], "C")
    atom_thr = Atom("N", "ALA", 'A', 20, [4.0, 9.0, 3.0], "N")
    atoms = [atom_one, atom_two, atom_thr]
    coords_test = atomcoords(atoms)
    coords_real = [
        1.0 8.0 4.0;
        4.0 4.0 9.0;
        5.0 2.0 3.0;
    ]
    @test coords_test == coords_real


    atom_one = Atom("CA", "ALA", 'A', 20, [1.0, 4.0, 5.0], "C")
    atom_two = Atom("C", "ALA", 'A', 20, [8.0, 4.0, 2.0], "C")
    atoms = [atom_one, atom_two]
    coords = [
        5.0 3.0;
        6.0 8.0;
        1.0 0.0;
    ]
    atoms_test = atomcoords(atoms, coords)
    @test atoms_test[1].atom_name == "CA"
    @test atoms_test[2].atom_name == "C"
    @test atoms_test[1].coords == coords[:,1]
    @test atoms_test[2].coords == coords[:,2]
    @test atoms[1].coords == [1.0, 4.0, 5.0]
    @test atoms[2].coords == [8.0, 4.0, 2.0]


    atom_one = Atom("CA", "ALA", 'A', 20, [1.0, 4.0, 5.0], "C")
    atom_two = Atom("C", "ALA", 'A', 20, [8.0, 4.0, 2.0], "C")
    atoms = [atom_one, atom_two]
    coords = [
        5.0 3.0;
        6.0 8.0;
        1.0 0.0;
    ]
    atomcoords!(atoms, coords)
    @test atoms[1].coords == coords[:,1]
    @test atoms[2].coords == coords[:,2]


    atom = Atom(mod_atom_info["atom_name"], mod_atom_info["res_name"], mod_atom_info["chain_id"], mod_atom_info["res_n"], [0.0, 0.0, 0.0], mod_atom_info["element"])
    @test ismodulator(atom)
    atom = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    @test !ismodulator(atom)
end
