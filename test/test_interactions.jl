# Tests for inters.jl


test_pdb_filepath = testfile("4AKE_H.pdb")
test_dssp_filepath = testfile("4AKE.dssp")


@testset "Interactions" begin
    atom_one = Atom("CA", "ALA", 'A', 20, [1.0, 0.0, 0.0], "C")
    atom_two = Atom("CA", "ALA", 'A', 21, [0.0, 1.0, 0.0], "C")
    atom_thr = Atom("CA", "ALA", 'A', 22, [0.0, 0.0, 1.0], "C")
    atoms = Atom[atom_one, atom_two, atom_thr]
    dists_test = distances(atoms)
    dists_real = [
        0.0   sqrt(2)   sqrt(2);
        sqrt(2)   0.0   sqrt(2);
        sqrt(2)   sqrt(2)   0.0;
    ]
    @test isapprox(dists_test, dists_real)


    bonded_pairs = BondedPair[
        BondedPair("ALA", "N", "CA"),
        BondedPair("ALA", "C", "O"),
        BondedPair("ALA", "C", "N")
    ]
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_thr = Atom("O", "ALA", 'A', 20, [0.0, 0.0, 0.0], "O")
    atoms = Atom[atom_one, atom_two, atom_thr]
    bonds_test = bonds(atoms, bonded_pairs)
    bonds_real = [
        false true false;
        true false false;
        false false false;
    ]
    @test bonds_test == bonds_real
    # Test if the amide bond is found correctly
    atom_one = Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("N", "ALA", 'A', 21, [0.0, 0.0, 0.0], "N")
    atom_thr = Atom("N", "ALA", 'A', 22, [0.0, 0.0, 0.0], "N")
    atoms = Atom[atom_one, atom_two, atom_thr]
    bonds_test = bonds(atoms, bonded_pairs)
    bonds_real = [
        false true false;
        true false false;
        false false false;
    ]
    @test bonds_test == bonds_real


    bonded_pairs = BondedPair[
        BondedPair("ALA", "N", "CA"),
        BondedPair("ALA", "CA", "C"),
    ]
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_thr = Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atoms = Atom[atom_one, atom_two, atom_thr]
    angles_test = angles(atoms, bonds(atoms, bonded_pairs))
    angles_real = [
        false false true;
        false false false;
        true false false;
    ]
    @test angles_test == angles_real
    # Test angles are detected across residues
    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_thr = Atom("N", "ALA", 'A', 21, [0.0, 0.0, 0.0], "N")
    atoms = Atom[atom_one, atom_two, atom_thr]
    angles_test = angles(atoms, bonds(atoms, bonded_pairs))
    angles_real = [
        false false true;
        false false false;
        true false false;
    ]
    @test angles_test == angles_real


    bonded_pairs = BondedPair[
        BondedPair("ALA", "N", "CA"),
        BondedPair("ALA", "CA", "C"),
    ]
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_three = Atom("C", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_four = Atom("N", "ALA", 'A', 21, [0.0, 0.0, 0.0], "N")
    atoms = [atom_one, atom_two, atom_three, atom_four]
    bs = bonds(atoms, bonded_pairs)
    divalents_test = divalents(atoms, bs, angles(atoms, bs))
    divalents_real = [
        false false false true;
        false false false false;
        false false false false;
        true false false false;
    ]
    @test divalents_test == divalents_real


    # It is assumed that atoms are different
    # Test if atoms in different residues return false
    atom_one = Atom("CD2", "HIS", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CG", "HIS", 'A', 30, [0.0, 0.0, 0.0], "C")
    @test isring(atom_one, atom_two) == false
    # Test if atoms on same residue but not both in ring return false
    atom_one = Atom("CD2", "HIS", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CA", "HIS", 'A', 20, [0.0, 0.0, 0.0], "C")
    @test isring(atom_one, atom_two) == false
    # Test if atoms on same residue and both in ring return true
    atom_one = Atom("CD2", "HIS", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CG", "HIS", 'A', 20, [0.0, 0.0, 0.0], "C")
    @test isring(atom_one, atom_two) == true


    # It is assumed that the atoms are a 1-4 pair
    # Test if atoms in a double bond as required return true
    atom_one = Atom("CB", "ASN", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("1HD1", "ASN", 'A', 20, [0.0, 0.0, 0.0], "H")
    @test isdbonefour(atom_one, atom_two) == true
    # Test if atoms not in a double bond as required return false
    atom_one = Atom("CG", "ASN", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("C", "ASN", 'A', 20, [0.0, 0.0, 0.0], "C")
    @test isdbonefour(atom_one, atom_two) == false


    # Assumes atoms are a 1-4 pair
    # Test that omega 1-4 bonds return true
    atom_one = Atom("O", "ALA", 'A', 20, [0.0, 0.0, 0.0], "O")
    atom_two = Atom("CA", "ALA", 'A', 21, [0.0, 0.0, 0.0], "CA")
    @test isomonefour(atom_one, atom_two) == true
    # Test that non omega 1-4 bonds return false
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("O", "ALA", 'A', 21, [0.0, 0.0, 0.0], "O")
    @test isomonefour(atom_one, atom_two) == false


    # Assumes atoms are a 1-4 pair
    # Test that phi/psi bonds return true
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("O", "ALA", 'A', 21, [0.0, 0.0, 0.0], "O")
    @test isphipsi(atom_one, atom_two) == true
    # Test that non phi/psi bonds return false
    atom_one = Atom("O", "ALA", 'A', 20, [0.0, 0.0, 0.0], "O")
    atom_two = Atom("CA", "ALA", 'A', 21, [0.0, 0.0, 0.0], "CA")
    @test isphipsi(atom_one, atom_two) == false


    # Assumes atoms are a 1-4 pair around the phi/psi angle
    # Test that atoms in the same helix return true
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("O", "ALA", 'A', 21, [0.0, 0.0, 0.0], "O")
    dssps = Dict{ASCIIString, Char}("20A"=> 'H', "21A"=> 'H')
    @test istightphipsi(atom_one, atom_two, dssps) == true
    # Test that atoms where one residue is proline return true
    atom_one = Atom("N", "PRO", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("O", "ALA", 'A', 21, [0.0, 0.0, 0.0], "O")
    dssps = Dict{ASCIIString, Char}("20A"=> ' ', "21A"=> ' ')
    @test istightphipsi(atom_one, atom_two, dssps) == true
    # Test that atoms not in a tight phi/psi bond return false
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("O", "ALA", 'A', 21, [0.0, 0.0, 0.0], "O")
    dssps = Dict{ASCIIString, Char}("20A"=> ' ', "21A"=> ' ')
    @test istightphipsi(atom_one, atom_two, dssps) == false


    # Assumes atoms are a 1-4 pair around the phi/psi angle
    # Test that atoms in a loop region return true
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("O", "ALA", 'A', 21, [0.0, 0.0, 0.0], "O")
    dssps = Dict{ASCIIString, Char}("20A"=> ' ', "21A"=> ' ')
    @test isloosephipsi(atom_one, atom_two, dssps) == true
    # Test that atoms where one residue is glycine return true
    atom_one = Atom("N", "GLY", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("O", "ALA", 'A', 21, [0.0, 0.0, 0.0], "O")
    dssps = Dict{ASCIIString, Char}("20A"=> 'H', "21A"=> 'H')
    @test isloosephipsi(atom_one, atom_two, dssps) == true
    # Test that atoms not in a loose phi/psi bond return false
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("O", "ALA", 'A', 21, [0.0, 0.0, 0.0], "O")
    dssps = Dict{ASCIIString, Char}("20A"=> 'H', "21A"=> 'H')
    @test isloosephipsi(atom_one, atom_two, dssps) == false


    # Test two backbone atoms in the same helix return true
    atoms = readpdb(test_pdb_filepath)
    dssps = readdssp(test_dssp_filepath, atoms)
    @test issecstr(atoms[145], atoms[178], atoms, dssps) == true
    # Test two backbone atoms in the same helix but more than 4 residues apart return false
    @test issecstr(atoms[145], atoms[210], atoms, dssps) == false
    # Test non-backbone atoms return false
    @test issecstr(atoms[145], atoms[181], atoms, dssps) == false
    # Test two backbone atoms not in the same secondary structural element return false
    @test issecstr(atoms[145], atoms[127], atoms, dssps) == false


    # Test if close groups of the same charge return false
    atom_one = Atom("CZ", "ARG", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("NZ", "LYS", 'A', 30, [3.0, 0.0, 0.0], "N")
    atoms = [atom_one, atom_two]
    @test issaltbridge(atom_one, atom_two, atoms, distances(atoms)) == false
    # Test if close groups of opposite charge return true
    atom_one = Atom("CZ", "ARG", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CG", "ASP", 'A', 30, [3.0, 0.0, 0.0], "C")
    atoms = [atom_one, atom_two]
    @test issaltbridge(atom_one, atom_two, atoms, distances(atoms)) == true
    # Test if groups of opposite charge with other close atoms return true
    atom_one = Atom("CZ", "ARG", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("CG", "ASP", 'A', 30, [5.0, 0.0, 0.0], "C")
    atom_three = Atom("NH1", "ARG", 'A', 20, [1.0, 0.0, 0.0], "N")
    atom_four = Atom("OD1", "ASP", 'A', 30, [4.0, 0.0, 0.0], "O")
    atoms = [atom_one, atom_two, atom_three, atom_four]
    @test issaltbridge(atom_one, atom_two, atoms, distances(atoms)) == true


    # Test donor and acceptor atoms far apart return false
    bonded_pairs = BondedPair[BondedPair("ALA", "N", "H")]
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("H", "ALA", 'A', 20, [1.0, 0.0, 0.0], "H")
    atom_three = Atom("O", "ALA", 'A', 30, [4.0, 0.0, 0.0], "O")
    atoms = [atom_one, atom_two, atom_three]
    @test ishbond(1, 3, atoms, distances(atoms), bonds(atoms, bonded_pairs)) == false
    # Test donor and acceptor atoms return true if H is at correct angle
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("H", "ALA", 'A', 20, [1.0, 0.0, 0.0], "H")
    atom_three = Atom("O", "ALA", 'A', 30, [3.0, 0.0, 0.0], "O")
    atoms = [atom_one, atom_two, atom_three]
    @test ishbond(1, 3, atoms, distances(atoms), bonds(atoms, bonded_pairs)) == true
    # Test donor and acceptor atoms return false if H is at wrong angle
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("H", "ALA", 'A', 20, [0.0, 1.0, 0.0], "H")
    atom_three = Atom("O", "ALA", 'A', 30, [3.0, 0.0, 0.0], "O")
    atoms = [atom_one, atom_two, atom_three]
    @test ishbond(1, 3, atoms, distances(atoms), bonds(atoms, bonded_pairs)) == false
    # Test acceptor and H atoms return true if they are close enough
    atom_one = Atom("N", "ALA", 'A', 20, [0.0, 0.0, 0.0], "N")
    atom_two = Atom("H", "ALA", 'A', 20, [1.0, 0.0, 0.0], "H")
    atom_three = Atom("O", "ALA", 'A', 30, [3.0, 0.0, 0.0], "O")
    atoms = [atom_one, atom_two, atom_three]
    @test ishbond(2, 3, atoms, distances(atoms), bonds(atoms, bonded_pairs)) == true


    # Test if hydrophobic atoms within distance return true
    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("H", "ALA", 'A', 30, [3.5, 0.0, 0.0], "H")
    @test ishydrophobic(atom_one, atom_two, 3.5, 1.0) == true
    # Test if hydrophobic atoms out of distance return false
    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("H", "ALA", 'A', 30, [4.5, 0.0, 0.0], "H")
    @test ishydrophobic(atom_one, atom_two, 4.5, 1.0) == false
    # Test if non-hydrophobic atoms in distance return false
    atom_one = Atom("CA", "ALA", 'A', 20, [0.0, 0.0, 0.0], "C")
    atom_two = Atom("O", "ALA", 'A', 30, [3.5, 0.0, 0.0], "O")
    @test ishydrophobic(atom_one, atom_two, 3.5, 1.0) == false


    atom_one = Atom("O", "ALA", 'A', 20, [0.0, 0.0, 0.0], "O")
    atom_two = Atom("C", "ALA", 'A', 20, [1.2, 0.0, 0.0], "C")
    atom_thr = Atom("H", "ALA", 'A', 30, [-1.0, -1.0, -1.0], "H")
    atoms = [atom_one, atom_two, atom_thr]
    dists = distances(atoms)
    inters = [
        0 1 12;
        1 0 15;
        12 15 0;
    ]
    bounds = Bounds(atoms, dists, inters, other_ratio=0.0, bound_weight=1.0)
    @test bounds.pres_inds == [2 1; 3 1]
    inter_types, tolerances = interactioninfo(1.0)
    lower_real = [
        dists[2, 1] - tolerances[1],
        dists[3, 1] - tolerances[12],
    ]
    upper_real = [
        dists[2, 1] + tolerances[1],
        dists[3, 1] + tolerances[12],
    ]
    @test isapprox(bounds.lower, lower_real)
    @test isapprox(bounds.upper, upper_real)


    atom_a_one = Atom("O", "ALA", 'A', 20, [0.0, 0.0, 0.0], "O")
    atom_a_two = Atom("C", "ALA", 'A', 20, [1.2, 0.0, 0.0], "C")
    atom_a_thr = Atom("H", "ALA", 'A', 30, [-1.0, -1.0, -1.0], "H")
    atoms_a = [atom_a_one, atom_a_two, atom_a_thr]
    atom_b_one = Atom("O", "ALA", 'A', 20, [0.0, 0.0, 0.0], "O")
    atom_b_two = Atom("C", "ALA", 'A', 20, [1.25, 0.0, 0.0], "C")
    atom_b_thr = Atom("H", "ALA", 'A', 30, [-1.1, -1.1, -1.1], "H")
    atoms_b = [atom_b_one, atom_b_two, atom_b_thr]
    dists_a = distances(atoms_a)
    dists_b = distances(atoms_b)
    inters_a = [
        0 1 12;
        1 0 15;
        12 15 0;
    ]
    inters_b = [
        0 1 10;
        1 0 15;
        10 15 0;
    ]
    bounds = Bounds(atoms_a, dists_a, dists_b, inters_a, inters_b, other_ratio=0.0, bound_weight=1.0)
    @test bounds.pres_inds == [2 1; 3 1]
    inter_types, tolerances = interactioninfo(1.0)
    lower_real = [
        min(dists_a[2, 1] - tolerances[1], dists_b[2, 1] - tolerances[1]),
        min(dists_a[3, 1] - tolerances[12], dists_b[3, 1] - tolerances[10]),
    ]
    upper_real = [
        max(dists_a[2, 1] + tolerances[1], dists_b[2, 1] + tolerances[1]),
        max(dists_a[3, 1] + tolerances[12], dists_b[3, 1] + tolerances[10]),
    ]
    @test isapprox(bounds.lower, lower_real)
    @test isapprox(bounds.upper, upper_real)
end
