# Read and write files


export
    checkfilepath,
    readdssp,
    readpdb,
    fixstring,
    spaceatomname,
    writepdb,
    readensemble,
    writeensemble,
    writeensemblescores,
    writeprojections,
    writepcview,
    writepcviews,
    writeintarray,
    writefloatarray,
    writestringarray,
    readpocketpoints,
    readligsite,
    readpdblines,
    writeclusterpoints


"Checks if an output filepath is valid, and throws an error if not."
function checkfilepath(out_filepath::AbstractString)
    directory = splitdir(out_filepath)[1]
    @assert directory == "" || isdir(directory) "Not a valid file location: \"$out_filepath\""
end


"""
Reads a DSSP file. If atoms are provided, creates additional '-' entries for residues not present in the DSSP file.
Returns a dictionary where key is residue number and chain, and value is secondary structure.
"""
function readdssp(dssp_filepath::AbstractString, atoms::Array{Atom,1})
    dssps = readdssp(dssp_filepath)
    res_pres = collect(keys(dssps))
    for atom in atoms
        res = "$(atom.res_n)$(atom.chain_id)"
        if !(res in res_pres)
            dssps[res] = '-'
            push!(res_pres, res)
        end
    end
    return dssps
end

function readdssp(dssp_filepath::AbstractString)
    expanded_path = expanduser(dssp_filepath)
    @assert isfile(expanded_path) "Not a valid filepath: \"$expanded_path\""
    dssps = Dict{String, Char}()
    reading = false
    open(expanded_path, "r") do dssp_file
        for line in eachline(dssp_file)
            if reading && length(line) >= 17 && line[14] != '!'
                res = "$(strip(line[7:10]))$(line[12])"
                sec_str = line[17]
                dssps[res] = sec_str
            elseif !reading && length(line) >= 12 && line[1:12] == "  #  RESIDUE"
                reading = true
            end
        end
    end
    if length(dssps) == 0
        println("Nothing read from DSSP file")
    else
        println("Read ", length(dssps), " residues from DSSP file")
    end
    return dssps
end


"""
Parse a PDB file.
Returns a list of atoms.
hetatm determines whether HETATM records are read in.
"""
function readpdb(in_filepath::AbstractString; hetatm::Bool=false)
    expanded_path = expanduser(in_filepath)
    @assert isfile(expanded_path) "Not a valid filepath: \"$expanded_path\""
    atoms = Atom[]
    open(expanded_path, "r") do in_file
        if hetatm
            line_starts = ["ATOM  ", "HETATM"]
        else
            line_starts = ["ATOM  "]
        end
        for line in eachline(in_file)
            if length(line) >= 54 && line[1:6] in line_starts
                # Atom name string is columns 13-16 with leading/trailing spaces removed
                atom_name = strip(line[13:16])
                # Residue name string is columns 18-20
                res_name = line[18:20]
                # Chain character is column 22
                chain_id = line[22]
                # Residue number integer is columns 23-26 converted to an integer
                res_n = Meta.parse(Int, line[23:26])
                # x, y and z coordinates are columns 31-38, 39-46 and 47-54 respectively converted to floats
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                # Element is determined from atom name, not from columns 77-78
                element = inferelement(atom_name)
                # Remove disorder
                if (line[17] == ' ' || line[17] == 'A')
                    push!(atoms, Atom(atom_name, res_name, chain_id, res_n, [x, y, z], element))
                end
            end
        end
    end
    println("Read ", length(atoms), " atoms from PDB file")
    return atoms
end


"""
Converts a string, character or number to a fixed length string.
Be careful - could e.g. cut 0.5 to "0".
Also does not round correctly.
"""
function fixstring(val_in, new_length::Integer; fill_char::Char=' ')
    string_out = string(val_in)
    string_out = string_out[1:min(length(string_out), new_length)]
    while length(string_out) < new_length
        string_out = "$fill_char$string_out"
    end
    return string_out
end


"""
Space an atom name correctly in the four-character limit.
The element is in position 2.
Returns a string of length 4.
"""
function spaceatomname(atom_name::AbstractString, element::AbstractString)
    chars = length(atom_name)
    @assert chars <= 4 "Atom name is greater than four characters"
    if element != "-"
        cent_ind = findfirst(atom_name, element[1])
    else
        cent_ind = 1
    end
    @assert cent_ind <= 2 "Atom name is too long to align correctly"
    @assert chars - cent_ind <= 3 "Atom name is too long to align correctly"
    if cent_ind == 1 && chars < 4
        out_string = " "
    else
        out_string = ""
    end
    out_string = "$out_string$atom_name"
    while length(out_string) < 4
        out_string = "$out_string "
    end
    return out_string
end


"""
Write a list of atoms as a PDB file.
The atoms are not reordered so it is assumed they are in a sensible order, e.g. the order they were read in.
Only ATOM records are written; no TER labels.
Could write a brief header for generated models.
"""
function writepdb(out_filepath::AbstractString, atoms::Array{Atom,1})
    expanded_path = expanduser(out_filepath)
    checkfilepath(expanded_path)
    open(expanded_path, "w") do out_file
        println(out_file, "REMARK   Generated by ExProSE on ", Libc.strftime(time()))
        for (atom_counter, atom) in enumerate(atoms)
            line = String[
                "ATOM  ",
                fixstring(atom_counter, 5),
                " ",
                spaceatomname(atom.atom_name, atom.element),
                " ",
                fixstring(atom.res_name, 3),
                " ",
                string(atom.chain_id),
                fixstring(atom.res_n, 4),
                "      ",
                fixstring(atom.coords[1], 6),
                "  ",
                fixstring(atom.coords[2], 6),
                "  ",
                fixstring(atom.coords[3], 6),
                "  1.00  0.00          ",
                fixstring(atom.element, 2),
                "  "
            ]
            @assert length(join(line)) == 80 "PDB line not correct length"
            println(out_file, line...)
        end
    end
end


"Read an ensemble of structures and return a ModelledEnsemble."
function readensemble(out_prefix::AbstractString, n_to_read::Integer)
    @assert n_to_read > 0 "n_to_read must be positive"
    atoms = readpdb("$out_prefix.1")
    strucs = ModelledStructure[]
    for i in 1:n_to_read
        curr_filepath = "$out_prefix.$i"
        @assert isfile(curr_filepath) "Not a valid filepath: \"$curr_filepath\""
        curr_atoms = readpdb(curr_filepath)
        # Placeholder score of -1.0 given to structures read in
        struc = ModelledStructure(-1.0, atomcoords(curr_atoms))
        push!(strucs, struc)
    end
    println("Read in protein ensemble")
    return ModelledEnsemble(atoms, strucs)
end


"Write an ensemble out as separate PDB files."
function writeensemble(out_prefix::AbstractString, ensemble::ModelledEnsemble)
    checkfilepath(out_prefix)
    atoms = ensemble.atoms
    strucs = ensemble.strucs
    @assert length(atoms) > 0 "No atoms associated with ensemble"
    @assert length(strucs) > 0 "The ensemble does not contain any structures"
    for (i, struc) in enumerate(strucs)
        out_filepath = "$out_prefix.$i"
        struc_atoms = atomcoords(atoms, struc.coords)
        writepdb(out_filepath, struc_atoms)
    end
    println("Wrote ensemble to PDB file(s) with prefix \"$out_prefix\"")
end


"Write out the scores associated with the structures in an ensemble."
function writeensemblescores(out_filepath::AbstractString, ensemble::ModelledEnsemble)
    expanded_path = expanduser(out_filepath)
    checkfilepath(expanded_path)
    strucs = ensemble.strucs
    @assert length(strucs) > 0 "The ensemble does not contain any structures"
    open(expanded_path, "w") do out_file
        for struc in strucs
            println(out_file, round(struc.score, digits=1))
        end
    end
    println("Wrote ensemble scores to file \"$expanded_path\"")
end


"""
Write out the projections of structures onto the eigenvectors of a PCA.
pcs_to_write is a list of the PCs to write, e.g. [1, 2].
All PCs written by default.
"""
function writeprojections(out_filepath::AbstractString, pcs::Array{Float64}, pcs_to_write::Array{Int,1}=collect(1:size(pcs, 1)))
    expanded_path = expanduser(out_filepath)
    checkfilepath(expanded_path)
    @assert size(pcs, 1) >= maximum(pcs_to_write) "Chose to write out a PC higher than the number of PCs present"
    open(expanded_path, "w") do out_file
        for i in 1:size(pcs, 2)
            vals_out = getindex(pcs[:,i], pcs_to_write)
            string_out = join([string(val) for val in vals_out], "\t")
            println(out_file, string_out)
        end
    end
    println("Projections written to file \"$expanded_path\"")
end


"Write a PyMol script to view the atom motions of a principal component."
function writepcview(out_filepath::AbstractString,
                    pca::PCA,
                    atoms::Array{Atom,1},
                    pc::Integer;
                    weighting::Real=50.0,
                    line_colour::AbstractString="magenta")
    expanded_path = expanduser(out_filepath)
    checkfilepath(expanded_path)
    @assert pc > 0 "Selected principal component must be positive"
    @assert pc <= size(pca.evecs, 2) "Selected principal component is too high"
    inds_to_use = calphaindices(atoms)
    @assert size(pca.evecs, 1) == 3 * length(inds_to_use) "PCA and atom list have different sizes"
    disps = pca.evecs[:,pc]
    open(expanded_path, "w") do out_file
        println(out_file, "set dash_gap, 0")
        println(out_file, "set dash_color, ", line_colour)
        for (i, j) in enumerate(inds_to_use)
            pc_coords = atoms[j].coords + disps[3i-2:3i] * weighting
            println(out_file, "pseudoatom pt", 2i-1, ", pos=[", round(atoms[j].coords[1], digits=2), ", ", round(atoms[j].coords[2], digits=2), ", ", round(atoms[j].coords[3], digits=2), "]")
            println(out_file, "pseudoatom pt", 2i, ", pos=[", round(pc_coords[1], digits=2), ", ", round(pc_coords[2], digits=2), ", ", round(pc_coords[3], digits=2), "]")
            println(out_file, "distance dist", i, ", /pt", 2i-1, ", /pt", 2i)
        end
        println(out_file, "hide everything, pt*")
        println(out_file, "hide labels")
    end
end


"Wrapper function for writepcview that produces visualisation scripts for the top `n_to_write` PCs."
function writepcviews(out_prefix::AbstractString,
                    pca::PCA,
                    atoms::Array{Atom,1},
                    n_to_write::Integer)
    expanded_prefix = expanduser(out_prefix)
    checkfilepath(expanded_prefix)
    if n_to_write <= size(pca.evecs, 2)
        n_to_use = n_to_write
    else
        n_to_use = size(pca.evecs, 2)
        println("n_to_write is larger than the number of eigenvectors - proceeding with all possible eigenvectors")
    end
    for i in 1:n_to_use
        out_filepath = "$expanded_prefix$i.pml"
        writepcview(out_filepath, pca, atoms, i)
    end
    println("PyMol script to view principal components written to file(s) \"$(expanded_prefix)x.pml\"")
end


"Write an array of integers to an output file."
function writeintarray(out_filepath::AbstractString, ints::Array{Int,1})
    expanded_path = expanduser(out_filepath)
    checkfilepath(expanded_path)
    open(expanded_path, "w") do out_file
        for i in ints
            println(out_file, i)
        end
    end
end


"Write an array of floats to an output file."
function writefloatarray(out_filepath::AbstractString, vals::Array{Float64,1}; dec_places::Integer=3)
    expanded_path = expanduser(out_filepath)
    checkfilepath(expanded_path)
    open(expanded_path, "w") do out_file
        for val in vals
            println(out_file, round(val, digits=dec_places))
        end
    end
end


"Write an array of strings to an output file."
function writestringarray(out_filepath::AbstractString, vals::Array{<:AbstractString,1})
    expanded_path = expanduser(out_filepath)
    checkfilepath(expanded_path)
    open(expanded_path, "w") do out_file
        for val in vals
            println(out_file, val)
        end
    end
end


"""
Read a custom pocket point file in PDB format where the residue number is the pocket number.
This file can be generated using a script from the LIGSITEcs pocket_r.pdb and pocket_all.pdb files.
Returns a dictionary where key is pocket number and value is point coords.
"""
function readpocketpoints(pdb_filepath::AbstractString)
    expanded_path = expanduser(pdb_filepath)
    @assert isfile(expanded_path) "Not a valid filepath: \"$expanded_path\""
    pock_points = Dict{Int, Array{Float64}}()
    counter = 0
    open(expanded_path, "r") do pdb_file
        for line in eachline(pdb_file)
            if startswith(line, "ATOM  ") || startswith(line, "HETATM")
                counter += 1
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                pocket_n = Meta.parse(Int, line[23:26])
                if haskey(pock_points, pocket_n)
                    pock_points[pocket_n] = hcat(pock_points[pocket_n], [x, y, z])
                else
                    pock_points[pocket_n] = [x, y, z]
                end
            end
        end
    end
    println("Read ", counter, " pocket points from pocket points PDB file")
    return pock_points
end


"""
Read a LIGSITEcs PDB output file.
Returns the pocket centre coordinate list and list of pocket volumes.
"""
function readligsite(pdb_filepath::AbstractString)
    expanded_path = expanduser(pdb_filepath)
    @assert isfile(expanded_path) "Not a valid filepath: \"$expanded_path\""
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    vols = Int[]
    counter = 0
    open(expanded_path, "r") do pdb_file
        for line in eachline(pdb_file)
            if startswith(line, "ATOM  ") || startswith(line, "HETATM")
                counter += 1
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                vol = Meta.parse(Int, line[23:26])
                push!(xs, x)
                push!(ys, y)
                push!(zs, z)
                push!(vols, vol)
            end
        end
    end
    coords = transpose(hcat(xs, ys, zs))
    println("Read ", counter, " pocket centres from LIGSITEcs output PDB file")
    return coords, vols
end


"Read the ATOM lines from a PDB file into an array."
function readpdblines(pdb_filepath::AbstractString)
    expanded_path = expanduser(pdb_filepath)
    @assert isfile(expanded_path) "Not a valid filepath: \"$expanded_path\""
    point_lines = String[]
    open(expanded_path, "r") do in_file
        for line in eachline(in_file)
            if startswith(line, "ATOM  ")
                push!(point_lines, chomp(line))
            end
        end
    end
    return point_lines
end


"Write output file with assignment in residue number column"
function writeclusterpoints(out_filepath::AbstractString,
                    point_lines::Array{<:AbstractString,1},
                    assignments::Array{Int,1})
    expanded_path = expanduser(out_filepath)
    checkfilepath(expanded_path)
    open(expanded_path, "w") do out_file
        println(out_file, "REMARK    Pocket number in residue number column")
        println(out_file, "REMARK    Pocket numbers found using cluster-ligsite")
        for (i, line) in enumerate(point_lines)
            print(out_file, line[1:22], lpad(string(assignments[i]), 4), line[27:end])
        end
    end
    println("Wrote pocket number assignments to file \"$expanded_path\"")
end
