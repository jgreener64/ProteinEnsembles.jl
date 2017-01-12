# Tests for perturb.jl


test_pocket_r_filepath = testfile("1DVM_pocket_r.pdb")
test_pocket_all_filepath = testfile("1DVM_pocket_all.pdb")


@testset "Perturb" begin
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("CA", "ALA", 'A', 20, [0.0, 1.0, 0.0], "C")
    atom_thr = Atom("C", "ALA", 'A', 20, [1.0, 1.0, 0.0], "C")
    atoms = [atom_one, atom_two, atom_thr]
    lower = [0.9, 0.9]
    upper = [1.1, 1.1]
    pres_inds = [
        2 1;
        3 2;
    ]
    constraints = Constraints(atoms, lower, upper, pres_inds)
    mod_coords = [
        0.5 1.5;
        -0.5 0.5;
        0.0 0.0;
    ]
    new_constraints = Constraints(atoms, constraints, mod_coords, intra_cutoff=1000.0,
        intra_tolerance=0.1, inter_cutoff=1.0, inter_tolerance=0.1)
    @test length(new_constraints.atoms) == 5
    @test new_constraints.atoms[4].res_name == mod_atom_info["res_name"]
    real_lower = [
        0.9,
        0.9,
        sqrt(2) - 0.1,
        sqrt(2) * 0.5 - 0.1,
        sqrt(2) * 0.5 - 0.1,
    ]
    real_upper = [
        1.1,
        1.1,
        sqrt(2) + 0.1,
        sqrt(2) * 0.5 + 0.1,
        sqrt(2) * 0.5 + 0.1,
    ]
    @test isapprox(new_constraints.lower, real_lower)
    @test isapprox(new_constraints.upper, real_upper)
    @test new_constraints.pres_inds == [
        2 1;
        3 2;
        5 4;
        4 1;
        5 3;
    ]


    pock_points = [
        1.0 2.0;
        3.0 4.0;
        5.0 6.0;
    ]
    mod_coords = repeatpocketpoints(pock_points)
    @test size(mod_coords) == (3,defaults["mod_n_points"])
    @test defaults["mod_n_points"] < 4 || mod_coords[:,1:4] == [
        1.0 2.0 1.0 2.0;
        3.0 4.0 3.0 4.0;
        5.0 6.0 5.0 6.0;
    ]
    pock_points = [1.0, 2.0, 3.0]
    mod_coords = repeatpocketpoints(pock_points)
    @test mod_coords == repeat([1.0; 2.0; 3.0], inner=[1,defaults["mod_n_points"]])
    pock_points = ones(3, 2*defaults["mod_n_points"])
    mod_coords = repeatpocketpoints(pock_points)
    @test mod_coords == ones(3, defaults["mod_n_points"])


    temp_filepath = tempname()
    clusterligsite(test_pocket_r_filepath, test_pocket_all_filepath, temp_filepath)
    rm(temp_filepath)
end
