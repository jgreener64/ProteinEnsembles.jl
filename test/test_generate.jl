# Tests for generate.jl


@testset "Generate" begin
    lower = [1.0]
    upper = [2.0]
    pres_inds = [1 3]
    constraints = Constraints(Atom[], lower, upper, pres_inds)
    coords = [
        0.0 0.0 2.5;
        0.0 1.0 0.0;
        0.0 0.0 0.0;
    ]
    @test errorscore(constraints, coords) == 0.25


    # Test L-amino acids are not inverted
    atom_ca = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_co = Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_n = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_cb = Atom("CB", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atoms = [atom_ca, atom_co, atom_n, atom_cb]
    coords_l = [
        0.0 1.0 -1.0 -1.0;
        0.0 -1.0 -1.0 1.0;
        0.0 1.0 -1.0 1.0;
    ]
    coords_fixed = [
        0.0 1.0 -1.0 -1.0;
        0.0 -1.0 -1.0 1.0;
        0.0 1.0 -1.0 1.0;
    ]
    fixchirality!(coords_l, atoms)
    @test coords_l == coords_fixed
    # Test D-amino acids are inverted
    coords_d = [
        0.0 -1.0 1.0 -1.0;
        0.0 -1.0 -1.0 1.0;
        0.0 -1.0 1.0 1.0;
    ]
    coords_fixed = [
        0.0 -1.0 1.0 -1.0;
        0.0 -1.0 -1.0 1.0;
        0.0 1.0 -1.0 -1.0;
    ]
    fixchirality!(coords_d, atoms)
    @test coords_d == coords_fixed


    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom(mod_atom_info["atom_name"], mod_atom_info["res_name"], mod_atom_info["chain_id"], mod_atom_info["res_n"], [0.0, 0.0, 0.0], mod_atom_info["element"])
    atom_three = Atom("CA", "ALA", 'A', 21, [0.0, 0.0, 0.0], "C")
    atom_four = Atom(mod_atom_info["atom_name"], mod_atom_info["res_name"], mod_atom_info["chain_id"], mod_atom_info["res_n"], [0.0, 0.0, 0.0], mod_atom_info["element"])
    atoms = [atom_one, atom_two, atom_three, atom_four]
    coords_one = [
        1.0 10.0 5.0 11.0;
        9.0 4.0 8.0 2.0;
        7.0 3.0 6.0 12.0;
    ]
    coords_two = [
        9.0 6.0 11.0 3.0;
        1.0 2.0 5.0 12.0;
        7.0 8.0 10.0 4.0;
    ]
    struc_one = ModelledStructure(10.0, coords_one)
    struc_two = ModelledStructure(20.0, coords_two)
    strucs = [struc_one, struc_two]
    ensemble = ModelledEnsemble(atoms, strucs)
    removemodulator!(ensemble)
    test_atoms = ensemble.atoms
    @test length(test_atoms) == 2
    @test test_atoms[1].atom_name == "CA"
    @test test_atoms[2].atom_name == "CA"
    test_strucs = ensemble.strucs
    @test strucs[1].score == 10.0
    real_coords_one = [
        1.0 5.0;
        9.0 8.0;
        7.0 6.0;
    ]
    real_coords_two = [
        9.0 11.0;
        1.0 5.0;
        7.0 10.0;
    ]
    @test strucs[1].coords == real_coords_one
    @test strucs[2].coords == real_coords_two


    struc_one = ModelledStructure(40.0, [0.0 0.0; 0.0 0.0; 0.0 0.0])
    struc_two = ModelledStructure(20.0, [0.0 0.0; 0.0 0.0; 0.0 0.0])
    struc_three = ModelledStructure(10.0, [0.0 0.0; 0.0 0.0; 0.0 0.0])
    struc_four = ModelledStructure(30.0, [0.0 0.0; 0.0 0.0; 0.0 0.0])
    strucs = [struc_one, struc_two, struc_three, struc_four]
    ensemble = ModelledEnsemble([], strucs)
    trimbyscore!(ensemble, 2)
    trim_strucs = ensemble.strucs
    @test length(trim_strucs) == 2
    @test trim_strucs[1].score == 10.0
    @test trim_strucs[2].score == 20.0
end
