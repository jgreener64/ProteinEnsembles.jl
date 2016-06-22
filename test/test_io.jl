# Tests for io.jl


test_pdb_filepath = testfile("1AKE.pdb")
test_dssp_filepath = testfile("4AKE.dssp")
test_pocket_points = testfile("1AKE_pocket_points.pdb")


@testset "IO" begin
    dssps = readdssp(test_dssp_filepath)
    @test dssps["21A"] == 'H'
    @test dssps["43A"] == ' '


    atoms = readpdb(test_pdb_filepath)
    @test isa(atoms, Array{Atom,1})
    @test length(atoms) == 3312
    atom = atoms[10]
    @test atom.atom_name == "CA"
    @test atom.res_n == 2
    @test atom.coords[2] == 49.969
    atoms = readpdb(test_pdb_filepath, hetatm=true)
    @test length(atoms) == 3804


    @test fixstring("Test string", 4) == "Test"
    @test fixstring('T', 4) == "   T"
    @test fixstring(0.123456, 5) == "0.123"
    @test fixstring(123, 5) == "  123"


    @test spaceatomname("CA", "C") == " CA "
    @test spaceatomname("CAB", "C") == " CAB"
    @test spaceatomname("1HB1", "H") == "1HB1"
    @test spaceatomname("CABA", "C") == "CABA"


    atoms = readpdb(test_pdb_filepath)
    # Temp file which is removed at the end
    temp_filepath = tempname()
    writepdb(temp_filepath, atoms)
    test_atoms = readpdb(temp_filepath)
    @test length(test_atoms) == 3312
    atom = test_atoms[10]
    @test atom.atom_name == "CA"
    @test atom.res_n == 2
    @test atom.coords[2] == 49.969
    rm(temp_filepath)


    pock_points = readpocketpoints(test_pocket_points)
    @test length(pock_points) == 18
    @test size(pock_points[2]) == (3,54)
    @test pock_points[6] == [
        37.08   37.08   37.08 ;
        42.827  42.827  43.827;
        31.49   32.49   31.49 ;
    ]
end
