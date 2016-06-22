# Deal with atoms


export
    inferelement,
    vdw_radius,
    hard_sphere_radius,
    areatomssame,
    findcommonatoms,
    atommap,
    atomid,
    readdssp,
    calphaindices,
    atomcoords,
    atomcoords!,
    mod_atom_info,
    ismodulator


"""
Determines the element of an atom from the atom name.
Returns the element as a string.
"""
function inferelement(atom_name::AbstractString)
    if contains(atom_name, "C")
        element = "C"
    elseif contains(atom_name, "N")
        element = "N"
    elseif contains(atom_name, "O")
        element = "O"
    elseif contains(atom_name, "S")
        element = "S"
    elseif contains(atom_name, "H")
        element = "H"
    elseif contains(atom_name, "P")
        element = "P"
    else
        println("Could not find element for atom with name ", atom_name)
        element = "-"
    end
    return element
end


"Dictionary of atomic VDW radii in Angstroms."
const vdw_radius = Dict{ASCIIString, Float64}(
    "C"=> 1.70,
    "N"=> 1.55,
    "O"=> 1.52,
    "S"=> 1.80,
    "H"=> 1.20,
    "P"=> 1.80,
    "-"=> 1.70, # Say unknown elements are like carbon
)


"Dictionary of atomic hard sphere radii in Angstroms."
const hard_sphere_radius = Dict{ASCIIString, Float64}(
    "C"=> 0.67,
    "N"=> 0.56,
    "O"=> 0.48,
    "S"=> 0.88,
    "H"=> 0.53,
    "P"=> 0.98,
    "-"=> 0.67, # Say unknown elements are like carbon
)


"""
Determines if two atoms are the same.
Atoms are the same if they have the same `atom_name`, `chain_id`, `res_n` and `element`.
`res_name` is not considered.
"""
function areatomssame(atom_one::Atom, atom_two::Atom)
    if atom_one.atom_name == atom_two.atom_name && atom_one.chain_id == atom_two.chain_id && atom_one.res_n == atom_two.res_n && atom_one.element == atom_two.element
        answer = true
    else
        answer = false
    end
    return answer
end


"""
Checks two atom arrays and finds common atoms.
Returns copies of input atom arrays with only common atoms remaining.
The returned arrays are therefore of equal length.
"""
function findcommonatoms(atoms_one::Array{Atom,1}, atoms_two::Array{Atom,1})
    common_atoms_one = Atom[]
    common_atoms_two = Atom[]
    for atom_one in atoms_one
        for atom_two in atoms_two
            if areatomssame(atom_one, atom_two)
                push!(common_atoms_one, atom_one)
                push!(common_atoms_two, atom_two)
                break
            end
        end
    end
    return common_atoms_one, common_atoms_two
end


"""
Finds the mapping of the atoms in the first list to those in the second.
Returns an array of length the same as the second list where the number
is the atom index in the first list.
Atoms missing from the first list are given 0.
"""
function atommap(atoms_one::Array{Atom,1}, atoms_two::Array{Atom,1})
    atom_map = Int[]
    for atom_two in atoms_two
        found = false
        for (i, atom_one) in enumerate(atoms_one)
            if areatomssame(atom_one, atom_two)
                push!(atom_map, i)
                found = true
                break
            end
        end
        if !found
            push!(atom_map, 0)
        end
    end
    @assert length(atom_map) == length(atoms_two)
    return atom_map
end


"A string describing an `Atom` in the form `atom_name`/`res_name`/`res_n`/`chain_id`."
atomid(atom::Atom) = "$(atom.atom_name)/$(atom.res_name)/$(atom.res_n)/$(atom.chain_id)"


"""
Reads a DSSP file, then creates additional '-' entries for residues not present in the DSSP file.
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


"Return as an array the indices in an array of atoms that correspond to C-alpha atoms."
function calphaindices(atoms::Array{Atom,1})
    ca_inds = Int[]
    for (i, atom) in enumerate(atoms)
        if atom.atom_name == "CA"
            push!(ca_inds, i)
        end
    end
    return ca_inds
end


"Return the list of coordinates from a list of atoms."
function atomcoords(atoms::Array{Atom,1})
    coords = zeros(3, length(atoms))
    for (i, atom) in enumerate(atoms)
        coords[:,i] = atom.coords
    end
    return coords
end


"""
Set the coordinates of a list of atoms from a list of coordinates.
Returns a new atoms list with updated coordinates.
"""
function atomcoords(atoms::Array{Atom,1}, coords::Array{Float64})
    new_atoms = deepcopy(atoms)
    atomcoords!(new_atoms, coords)
    return new_atoms
end


"Set the coordinates of a list of atoms from a list of coordinates."
function atomcoords!(atoms::Array{Atom,1}, coords::Array{Float64})
    n_atoms = length(atoms)
    @assert n_atoms == size(coords, 2) "The number of atoms is not the same as the number of coordinates"
    for (i, atom) in enumerate(atoms)
        atom.coords = coords[:,i]
    end
end


"Dictionary of atomic information for fake modulator atoms."
const mod_atom_info = Dict{ASCIIString, Any}(
    "atom_name"=> "X",
    "res_name"=> "MOD",
    "chain_id"=> 'X',
    "res_n"=> 0,
    "element"=> "X",
)


"Determines whether an `Atom` is a fake modulator atom."
function ismodulator(atom::Atom)
    mod_atom = Atom(mod_atom_info["atom_name"], mod_atom_info["res_name"], mod_atom_info["chain_id"], mod_atom_info["res_n"], [0.0, 0.0, 0.0], mod_atom_info["element"])
    if areatomssame(atom, mod_atom)
        answer = true
    else
        answer = false
    end
    return answer
end
