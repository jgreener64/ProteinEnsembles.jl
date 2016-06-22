# Tests for perturb.jl


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
    bounds = Bounds(atoms, lower, upper, pres_inds)
    mod_coords = [
        0.5 1.5;
        -0.5 0.5;
        0.0 0.0;
    ]
    new_bounds = Bounds(atoms, bounds, mod_coords, intra_cutoff=1000.0,
        intra_freedom=0.1, inter_cutoff=1.0, inter_freedom=0.1)
    @test length(new_bounds.atoms) == 5
    @test new_bounds.atoms[4].res_name == mod_atom_info["res_name"]
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
    @test isapprox(new_bounds.lower, real_lower)
    @test isapprox(new_bounds.upper, real_upper)
    @test new_bounds.pres_inds == [
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
    @test size(mod_coords) == (3,defaults["mod_points"])
    @test defaults["mod_points"] < 4 || mod_coords[:,1:4] == [
        1.0 2.0 1.0 2.0;
        3.0 4.0 3.0 4.0;
        5.0 6.0 5.0 6.0;
    ]
    pock_points = [1.0, 2.0, 3.0]
    mod_coords = repeatpocketpoints(pock_points)
    @test mod_coords == repeat([1.0; 2.0; 3.0], inner=[1,defaults["mod_points"]])
    pock_points = ones(3, 2*defaults["mod_points"])
    mod_coords = repeatpocketpoints(pock_points)
    @test mod_coords == ones(3, defaults["mod_points"])
end
