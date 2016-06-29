# Determine the interactions between atom pairs in a PDB file


export
    distances,
    bonds,
    angles,
    divalents,
    isring,
    isdbonefour,
    isomonefour,
    isphipsi,
    istightphipsi,
    isloosephipsi,
    issecstr,
    issaltbridge,
    ishbond,
    ishydrophobic,
    interactioninfo,
    findinteractions,
    interactions


"""
Calculate the distance matrix for a list of atoms.
Returns a 2D array of distances.
"""
function distances(atoms::Array{Atom,1})
    n_atoms = length(atoms)
    @assert n_atoms > 0 "No atoms in atom list"
    dists = zeros(n_atoms, n_atoms)
    for i in 1:n_atoms
        for j in 1:i-1
            dist = norm(atoms[i].coords - atoms[j].coords)
            dists[i,j] = dist
            dists[j,i] = dist
        end
    end
    return dists
end


"""
Calculate the bonding matrix for a list of atoms.
Returns a 2D boolean array of bonding.
"""
function bonds(atoms::Array{Atom,1}, bonded_pairs::Array{BondedPair,1})
    n_atoms = length(atoms)
    @assert n_atoms > 0 "No atoms in atom list"
    n_cov_dics = length(bonded_pairs)
    @assert n_cov_dics > 0 "No bonds in bonding dictionary"
    bs = falses(n_atoms, n_atoms)
    n_bonds = 0
    for i in 1:n_atoms
        for j in 1:i-1
            bonded = false
            # Check if the atoms are in the same residue
            if atoms[i].res_n == atoms[j].res_n && atoms[i].chain_id == atoms[j].chain_id
                # Check if the atoms appear as a bond in the bonding dictionary
                for k in 1:n_cov_dics
                    if atoms[i].res_name == bonded_pairs[k].residue
                        if atoms[i].atom_name == bonded_pairs[k].atom_one && atoms[j].atom_name == bonded_pairs[k].atom_two
                            bonded = true
                            break
                        elseif atoms[i].atom_name == bonded_pairs[k].atom_two && atoms[j].atom_name == bonded_pairs[k].atom_one
                            bonded = true
                            break
                        end
                    end
                end
            # Check if the atoms are in adjacent residues and part of the amide bond
            elseif atoms[i].chain_id == atoms[j].chain_id && abs(atoms[i].res_n - atoms[j].res_n) == 1
                if atoms[i].res_n > atoms[j].res_n && atoms[i].atom_name == "N" && atoms[j].atom_name == "C"
                    bonded = true
                elseif atoms[j].res_n > atoms[i].res_n && atoms[i].atom_name == "C" && atoms[j].atom_name == "N"
                    bonded = true
                end
            end
            bs[i,j] = bonded
            bs[j,i] = bonded
            if bonded
                n_bonds += 1
            end
        end
    end
    return bs
end


"""
Calculate the bond angle matrix for a list of atoms.
Returns a 2D boolean array of whether bond angles are present.
"""
function angles(atoms::Array{Atom,1}, bonds::BitArray)
    n_atoms = length(atoms)
    @assert n_atoms > 0 "No atoms in atom list"
    @assert n_atoms == size(bonds, 1) "Number of atoms in atom list and bonding matrix are not the same"
    angles = falses(n_atoms, n_atoms)
    n_angles = 0
    for i in 1:n_atoms
        for j in 1:i-1
            # Atoms have to be on the same or adjacent residues
            if atoms[i].chain_id == atoms[j].chain_id && abs(atoms[i].res_n - atoms[j].res_n) <= 1
                # Covalently-bonded atoms cannot be part of a 1-3 system
                if !bonds[i, j]
                    for k in 1:n_atoms
                        if bonds[i, k] && bonds[k, j]
                            angles[i, j] = true
                            angles[j, i] = true
                            n_angles += 1
                            break
                        end
                    end
                end
            end
        end
    end
    return angles
end


"""
Calculate the divalent matrix for a list of atoms.
Returns a 2D boolean array of whether divalents are present.
"""
function divalents(atoms::Array{Atom,1}, bonds::BitArray, angles::BitArray)
    n_atoms = length(atoms)
    @assert n_atoms > 0 "No atoms in atom list"
    @assert n_atoms == size(bonds, 1) "Number of atoms in atom list and bonding matrix are not the same"
    @assert n_atoms == size(angles, 1) "Number of atoms in atom list and angle matrix are not the same"
    divalents = falses(n_atoms, n_atoms)
    n_divalents = 0
    for i in 1:n_atoms
        for j in 1:i-1
            # Atoms have to be on the same or adjacent residues
            if atoms[i].chain_id == atoms[j].chain_id && abs(atoms[i].res_n - atoms[j].res_n) <= 1
                # Being part of a covalent or 1-3 system prevents a pair from being part of a 1-4 system
                if !bonds[i, j] && !angles[i, j]
                    for k in 1:n_atoms
                        if angles[i, k] && bonds[k, j]
                            divalents[i, j] = true
                            divalents[j, i] = true
                            n_divalents += 1
                            break
                        end
                    end
                end
            end
        end
    end
    return divalents
end


"""
Determines if two atoms are part of the same ring system.
"""
function isring(atom_one::Atom, atom_two::Atom)
    answer = false
    # See if the atoms are in the same residue
    if atom_one.res_n == atom_two.res_n && atom_one.chain_id == atom_two.chain_id
        if atom_one.res_name == "HIS"
            ring_res = ["CD2", "CG", "CE1", "ND1", "NE2", "1HD1", "1HD2", "1HE1", "1HE2"]
        elseif atom_one.res_name == "PHE"
            ring_res = ["CD2", "CG", "CZ", "CD1", "CE1", "CE2", "HD1", "HD2", "HE1", "HE2", "HZ"]
        elseif atom_one.res_name == "TYR"
            ring_res = ["CD2", "CG", "CZ", "CD1", "CE1", "CE2", "HD1", "HD2", "HE1", "HE2"]
        elseif atom_one.res_name == "TRP"
            ring_res = ["CZ2", "CZ3", "CG", "CH2", "CE3", "CE2", "CD2", "NE1", "CD1", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2"]
        else
            ring_res = []
        end
        if atom_one.atom_name in ring_res && atom_two.atom_name in ring_res
            answer = true
        end
    end
    return answer
end


"""
Determines if two atoms are a 1-4 pair around a restricted side chain double bond.
Assumes atoms are a 1-4 pair.
"""
function isdbonefour(atom_one::Atom, atom_two::Atom)
    answer = false
    # See if the atoms are in the same residue
    if atom_one.res_n == atom_two.res_n && atom_one.chain_id == atom_two.chain_id
        if atom_one.res_name == "ASN"
            rest_res = ["CB", "OD1", "1HD1", "2HD1"]
        elseif atom_one.res_name == "GLN"
            rest_res = ["CG", "OE1", "1HE2", "2HE2"]
        elseif atom_one.res_name == "ARG"
            rest_res = ["CD", "NE", "NH1", "NH2", "HE", "1HH1", "2HH1", "1HH2", "2HH2"]
        else
            rest_res = []
        end
        if atom_one.atom_name in rest_res && atom_two.atom_name in rest_res
            answer = true
        end
    end
    return answer
end


"""
Determines if two atoms are a 1-4 pair around the omega angle.
Assumes atoms are a 1-4 pair.
"""
function isomonefour(atom_one::Atom, atom_two::Atom)
    om_res = ["CA", "O", "H"]
    if atom_one.atom_name in om_res && atom_two.atom_name in om_res
        answer = true
    else
        answer = false
    end
    return answer
end


"""
Determines if two atoms are a 1-4 pair around the phi/psi angle.
Assumes atoms are a 1-4 pair.
Glycine H is not polar so should not be present.
"""
function isphipsi(atom_one::Atom, atom_two::Atom)
    pp_res = ["N", "O", "CB", "C", "H"]
    if atom_one.atom_name in pp_res && atom_two.atom_name in pp_res
        answer = true
    else
        answer = false
    end
    return answer
end


"""
Determines if two atoms in a 1-4 pair around the phi/psi angle are either in the same helix/strand or one of the residues is proline.
Assumes atoms are a 1-4 pair around the phi/psi angle.
"""
function istightphipsi(atom_one::Atom, atom_two::Atom, dssps::Dict{ASCIIString, Char})
    answer = false
    if atom_one.res_name == "PRO" || atom_two.res_name == "PRO"
        answer = true
    else
        dssp_one = dssps["$(atom_one.res_n)$(atom_one.chain_id)"]
        dssp_two = dssps["$(atom_two.res_n)$(atom_two.chain_id)"]
        if dssp_one == dssp_two && (dssp_one == 'H' || dssp_one == 'E')
            answer = true
        end
    end
    return answer
end


"""
Determines if two atoms in a 1-4 pair around the phi/psi angle are either in a loop region or one of the residues is glycine.
Assumes atoms are a 1-4 pair around the phi/psi angle.
"""
function isloosephipsi(atom_one::Atom, atom_two::Atom, dssps::Dict{ASCIIString, Char})
    answer = false
    if atom_one.res_name == "GLY" || atom_two.res_name == "GLY"
        answer = true
    else
        dssp_one = dssps["$(atom_one.res_n)$(atom_one.chain_id)"]
        dssp_two = dssps["$(atom_two.res_n)$(atom_two.chain_id)"]
        if dssp_one == dssp_two && dssp_one == ' '
            answer = true
        end
    end
    return answer
end


"""
Determines if two atoms are backbone atoms that are part of the same helix/strand and not more than 4 residues apart.
Does not include 3-10 and pi helices.
"""
function issecstr(atom_one::Atom, atom_two::Atom, atoms::Array{Atom,1}, dssps::Dict{ASCIIString, Char})
    answer = false
    if abs(atom_one.res_n - atom_two.res_n) <= 4 && atom_one.chain_id == atom_two.chain_id
        # See if the atoms have the same secondary structure that is a helix/strand
        dssp_one = dssps["$(atom_one.res_n)$(atom_one.chain_id)"]
        dssp_two = dssps["$(atom_two.res_n)$(atom_two.chain_id)"]
        if dssp_one == dssp_two && (dssp_one == 'H' || dssp_one == 'E')
            # See if the atoms are part of the backbone - currently not counting amide H
            bb_atoms = ["CA", "C", "O", "N"]
            if atom_one.atom_name in bb_atoms && atom_two.atom_name in bb_atoms
                answer = true
            end
            # Check residues between are in the same helix/strand
            for i in 1:length(atoms)
                if atoms[i].res_n > min(atom_one.res_n, atom_two.res_n) && atoms[i].res_n < max(atom_one.res_n, atom_two.res_n) && dssps["$(atoms[i].res_n)$(atoms[i].chain_id)"] != dssp_one && atoms[i].chain_id == atom_one.chain_id
                    answer = false
                    break
                end
            end
        end
    end
    return answer
end


"""
Determines if two atoms are part of oppositely charged groups that are within 4 Angstroms.
His not included.
"""
function issaltbridge(atom_one::Atom, atom_two::Atom, atoms::Array{Atom,1}, dists::Array{Float64})
    answer = false
    # See if the atoms are part of oppositely-charged residues
    if (atom_one.res_name == "ARG" || atom_one.res_name == "LYS") && (atom_two.res_name == "ASP" || atom_two.res_name == "GLU")
        atom_pos = atom_one
        atom_neg = atom_two
    elseif (atom_one.res_name == "ASP" || atom_one.res_name == "GLU") && (atom_two.res_name == "ARG" || atom_two.res_name == "LYS")
        atom_neg = atom_one
        atom_pos = atom_two
    else
        return false
    end
    # See if the atoms are found in the salt bridge parts of their residues
    if atom_pos.res_name == "ARG"
        opts_pos = ["CZ", "NH1", "NH2", "1HH1", "2HH1", "1HH2", "2HH2"]
    elseif atom_pos.res_name == "LYS"
        opts_pos = ["NZ"]
    end
    if atom_neg.res_name == "ASP"
        opts_neg = ["CG", "OD1", "OD2"]
    elseif atom_neg.res_name == "GLU"
        opts_neg = ["CD", "OE1", "OE2"]
    end
    if atom_pos.atom_name in opts_pos && atom_neg.atom_name in opts_neg
        # See if the closest interatomic distance between the positive and negative parts is close enough
        n_atoms = length(atoms)
        ind_pos = Int[]
        ind_neg = Int[]
        for i in 1:n_atoms
            if atom_pos.res_n == atoms[i].res_n && atom_pos.chain_id == atoms[i].chain_id && atoms[i].atom_name in opts_pos
                push!(ind_pos, i)
            elseif atom_neg.res_n == atoms[i].res_n && atom_neg.chain_id == atoms[i].chain_id && atoms[i].atom_name in opts_neg
                push!(ind_neg, i)
            end
        end
        for i in ind_pos
            for j in ind_neg
                if dists[i,j] < 4.0
                    answer = true
                    break
                end
            end
        end
    end
    return answer
end


"""
Determines if two atoms are part of a hydrogen bond where the donor-acceptor distance is <= 3.5 Angstroms, the
hydrogen-acceptor distance is <= 2.5 Angstroms and the donor-hydrogen-acceptor angle is minimally 90 degrees.
Won't pick up the donor-H pairing, even though it is part of the H bond, as this pairing is covalent.
Sulfur not considered.
"""
function ishbond(atom_id_one::Integer, atom_id_two::Integer, atoms::Array{Atom,1}, dists::Array{Float64}, bonds::BitArray)
    answer = false
    # Apply crude distance limit and check that atoms are not on the same residue
    if dists[atom_id_one, atom_id_two] <= 3.5 && (atoms[atom_id_one].res_n != atoms[atom_id_two].res_n || atoms[atom_id_one].chain_id != atoms[atom_id_two].chain_id)
        pol_res = ["N", "O"]
        n_atoms = length(atoms)
        atom_h_ids = Int[]
        # Crude way of re-checking distance later
        dist = 100
        # Case where the two atoms are the acceptor and donor
        if atoms[atom_id_one].element in pol_res && atoms[atom_id_two].element in pol_res
            atom_id_a = atom_id_one
            atom_id_b = atom_id_two
            # Find all Hs bonded to either atom that are within 2.5 Angstroms of the other atom
            for i in 1:n_atoms
                if atoms[i].element == "H" && ((bonds[i,atom_id_one] && dists[i,atom_id_two] <= 2.5) || (bonds[i,atom_id_two] && dists[i,atom_id_one] <= 2.5))
                    push!(atom_h_ids, i)
                end
            end
            # This test was already passed above
            dist = 0
        # Case where atom one is acceptor and atom two is H
        elseif dists[atom_id_one,atom_id_two] <= 2.5 && atoms[atom_id_one].element in pol_res && atoms[atom_id_two].element == "H"
            push!(atom_h_ids, atom_id_two)
            atom_id_a = atom_id_one
            for i in 1:n_atoms
                if bonds[i,atom_id_two]
                    atom_id_b = i
                    dist = dists[i,atom_id_one]
                    break
                end
            end
        # Case where atom two is acceptor and atom one is H
        elseif dists[atom_id_one,atom_id_two] <= 2.5 && atoms[atom_id_two].element in pol_res && atoms[atom_id_one].element == "H"
            push!(atom_h_ids, atom_id_one)
            atom_id_a = atom_id_two
            for i in 1:n_atoms
                if bonds[i,atom_id_one]
                    atom_id_b = i
                    dist = dists[i,atom_id_two]
                    break
                end
            end
        else
            return false
        end
        if dist <= 3.5
            for atom_id_h in atom_h_ids
                vec_a = atoms[atom_id_h].coords - atoms[atom_id_a].coords
                vec_b = atoms[atom_id_h].coords - atoms[atom_id_b].coords
                if dot(vec_a, vec_b) <= 0
                    answer = true
                    break
                end
            end
        end
    end
    return answer
end


"""
Determines if two atoms are closer than the sum of the VDW radii plus extra_dist Angstroms.
Currently counting C and H atoms only.
Applies a crude 5 A distance cutoff for speed, so extra_dist values above ~1.5 A will have no effect.
"""
function ishydrophobic(atom_one::Atom, atom_two::Atom, distance::Real, extra_dist::Real)
    if distance < 5 && distance < (vdw_radius[atom_one.element] + vdw_radius[atom_two.element] + extra_dist)
        hydr_list = ["C", "H"]
        if atom_one.element in hydr_list && atom_two.element in hydr_list
            answer = true
        else
            answer = false
        end
    else
        answer = false
    end
    return answer
end


"""
Returns a list with interaction types and a list with corresponding distance tolerances.
Zero index represents same atom.
`tolerance_weight` weights all tolerances.
"""
function interactioninfo(tolerance_weight::Real=defaults["tolerance_weight"])
    @assert tolerance_weight >= 0.0 "tolerance_weight cannot be negative"
    n_inter_types = 15
    inter_types = repeat([""], inner=[n_inter_types])
    tolerances = zeros(n_inter_types)
    inter_types[1] = "Covalently bonded"
    tolerances[1] = 0.02
    inter_types[2] = "1-3 system"
    tolerances[2] = 0.05
    inter_types[3] = "Same ring system"
    tolerances[3] = 0.1
    inter_types[4] = "Side-chain restricted 1-4"
    tolerances[4] = 0.1
    inter_types[5] = "Omega 1-4"
    tolerances[5] = 0.1
    inter_types[6] = "Tight phi/psi 1-4"
    tolerances[6] = 0.2
    inter_types[7] = "Loose phi/psi 1-4"
    tolerances[7] = 0.4
    inter_types[8] = "Other phi/psi 1-4"
    tolerances[8] = 0.3
    inter_types[9] = "Other 1-4"
    tolerances[9] = 0.4
    inter_types[10] = "Secondary structure"
    tolerances[10] = 0.5
    inter_types[11] = "Salt bridge"
    tolerances[11] = 0.75
    inter_types[12] = "H bond"
    tolerances[12] = 0.5
    inter_types[13] = "Tight hydrophobic"
    tolerances[13] = 0.5
    inter_types[14] = "Loose hydrophobic"
    tolerances[14] = 1.0
    inter_types[15] = "Other pairs"
    tolerances[15] = 5.0
    tolerances *= tolerance_weight
    return inter_types, tolerances
end


"""
Find the interactions between atoms in a list of atoms.
Returns 2D array of interactions by integer label.
Interactions the same as De Groot et. al., Structure, 1997, but cis/trans limits not implemented.
"""
function findinteractions(atoms::Array{Atom,1},
                    dists::Array{Float64},
                    bonds::BitArray,
                    angles::BitArray,
                    divalents::BitArray,
                    dssps::Dict{ASCIIString, Char})
    n_atoms = length(atoms)
    @assert n_atoms > 0 "No atoms in atom list"
    @assert n_atoms == size(dists, 1) "Number of atoms in atom list and distance matrix are not the same"
    @assert n_atoms == size(bonds, 1) "Number of atoms in atom list and bonding matrix are not the same"
    @assert n_atoms == size(angles, 1) "Number of atoms in atom list and angle matrix are not the same"
    @assert n_atoms == size(divalents, 1) "Number of atoms in atom list and divalent matrix are not the same"
    @assert length(dssps) > 0 "DSSP dictionary is empty"
    inters = zeros(Int, n_atoms, n_atoms)
    inter_types, tolerances = interactioninfo()
    inter_counter = zeros(Int, length(inter_types))
    for i in 1:n_atoms
        for j in 1:i-1
            if bonds[i,j]
                n_to_set = 1
            elseif angles[i,j]
                n_to_set = 2
            elseif isring(atoms[i], atoms[j])
                n_to_set = 3
            elseif divalents[i,j]
                if isdbonefour(atoms[i], atoms[j])
                    n_to_set = 4
                elseif isomonefour(atoms[i], atoms[j])
                    n_to_set = 5
                elseif isphipsi(atoms[i], atoms[j])
                    if istightphipsi(atoms[i], atoms[j], dssps)
                        n_to_set = 6
                    elseif isloosephipsi(atoms[i], atoms[j], dssps)
                        n_to_set = 7
                    else
                        n_to_set = 8
                    end
                else
                    n_to_set = 9
                end
            elseif issecstr(atoms[i], atoms[j], atoms, dssps)
                n_to_set = 10
            elseif issaltbridge(atoms[i], atoms[j], atoms, dists)
                n_to_set = 11
            elseif ishbond(i, j, atoms, dists, bonds)
                n_to_set = 12
            elseif ishydrophobic(atoms[i], atoms[j], dists[i, j], 0.5)
                n_to_set = 13
            elseif ishydrophobic(atoms[i], atoms[j], dists[i, j], 1.0)
                n_to_set = 14
            else
                n_to_set = 15
            end
            inters[i,j] = n_to_set
            inters[j,i] = n_to_set
            inter_counter[n_to_set] += 1
        end
    end
    println("Found interactions:")
    max_length = maximum(map(length, inter_types))
    for (k, inter_type) in enumerate(inter_types)
        println("  ", rpad(inter_type, max_length), " - ", inter_counter[k])
    end
    return inters
end


"""
Finds distance constraints from a set of interactions and distances.
Returns a `Constraints` object.
All explicit interactions calculated - other interactions written with probability such that `other_ratio` times n_atoms other interactions are written.
Sparse interactions scale as N but other interactions scale as N^2, hence the use of a ratio
rather than a direct probability to keep the sparse/other ratio roughly constant.
The hard sphere radii are multiplied by 0.8 to reduce the number of erroneous constraints.
"""
function Constraints(atoms::Array{Atom,1},
                    dists::Array{Float64},
                    inters::Array{Int};
                    other_ratio::Real=defaults["other_ratio"],
                    tolerance_weight::Real=defaults["tolerance_weight"])
    n_atoms = length(atoms)
    @assert n_atoms > 0 "No atoms in atom list"
    @assert n_atoms == size(dists, 1) "Number of atoms in atom list and distance matrix are not the same"
    @assert n_atoms == size(inters, 1) "Number of atoms in atom list and interaction matrix are not the same"
    @assert 0.0 <= other_ratio "other_ratio must be positive or zero"
    constraints_lower = Float64[]
    constraints_upper = Float64[]
    inter_types, tolerances = interactioninfo(tolerance_weight)
    present_i = Int[]
    present_j = Int[]
    # Calculate the probability of accepting other interactions based on the desired ratio of other interactions to n_atoms
    other_prob = min((2 * other_ratio) / (n_atoms - 1), 1.0)
    for i in 1:n_atoms
        for j in 1:i-1
            # Calculate explicit interactions and other interactions with a probability other_prob
            if inters[i, j] != 15 || rand() < other_prob
                constraint_upper = dists[i,j] + tolerances[inters[i,j]]
                # Factor of 0.8 used here
                radius_weighting = 0.8
                sum_radii = (hard_sphere_radius[atoms[i].element] + hard_sphere_radius[atoms[j].element]) * radius_weighting
                constraint_lower = max(dists[i,j] - tolerances[inters[i,j]], sum_radii)
                if constraint_lower > constraint_upper
                    println(
                        "Erroneous constraint:\t",
                        atomid(atoms[i]),
                        "\t",
                        atomid(atoms[j]),
                        "\t",
                        round(constraint_lower, 3),
                        "\t",
                        round(constraint_upper, 3)
                    )
                    constraint_lower = constraint_upper
                end
                push!(constraints_lower, constraint_lower)
                push!(constraints_upper, constraint_upper)
                push!(present_i, i)
                push!(present_j, j)
            end
        end
    end
    println("Found distance constraints from single structure")
    return Constraints(atoms, constraints_lower, constraints_upper, hcat(present_i, present_j))
end


"""
Finds distance constraints from two sets of interactions and distances.
Returns a `Constraints` object.
Input atoms must correspond directly to each other.
"""
function Constraints(atoms::Array{Atom,1},
                    dists_one::Array{Float64},
                    dists_two::Array{Float64},
                    inters_one::Array{Int},
                    inters_two::Array{Int};
                    other_ratio::Real=defaults["other_ratio"],
                    tolerance_weight::Real=defaults["tolerance_weight"])
    n_atoms = length(atoms)
    @assert n_atoms > 0 "No atoms in atom list"
    @assert n_atoms == size(dists_one, 1) "Number of atoms in atom list and distance matrix one are not the same"
    @assert n_atoms == size(dists_two, 1) "Number of atoms in atom list and distance matrix two are not the same"
    @assert n_atoms == size(inters_one, 1) "Number of atoms in atom list and interaction matrix one are not the same"
    @assert n_atoms == size(inters_two, 1) "Number of atoms in atom list and interaction matrix two are not the same"
    @assert 0.0 <= other_ratio "other_ratio must be positive or zero"
    constraints_lower = Float64[]
    constraints_upper = Float64[]
    inter_types, tolerances = interactioninfo(tolerance_weight)
    present_i = Int[]
    present_j = Int[]
    # Calculate the probability of accepting other interactions based on the desired ratio of other interactions to n_atoms
    other_prob = min((2 * other_ratio) / (n_atoms - 1), 1.0)
    for i in 1:n_atoms
        for j in 1:i-1
            # Calculate explicit interactions and other interactions with a probability other_prob
            if (inters_one[i, j] != 15 && inters_two[i, j] != 15) || rand() < other_prob
                constraint_upper_one = dists_one[i,j] + tolerances[inters_one[i,j]]
                constraint_upper_two = dists_two[i,j] + tolerances[inters_two[i,j]]
                constraint_upper = max(constraint_upper_one, constraint_upper_two)
                # Factor of 0.8 used here
                radius_weighting = 0.8
                sum_radii = (hard_sphere_radius[atoms[i].element] + hard_sphere_radius[atoms[j].element]) * radius_weighting
                constraint_lower_one = dists_one[i,j] - tolerances[inters_one[i,j]]
                constraint_lower_two = dists_two[i,j] - tolerances[inters_two[i,j]]
                constraint_lower = max(min(constraint_lower_one, constraint_lower_two), sum_radii)
                if constraint_lower > constraint_upper
                    println(
                        "Erroneous constraint:\t",
                        atomid(atoms[i]),
                        "\t",
                        atomid(atoms[j]),
                        "\t",
                        round(constraint_lower, 3),
                        "\t",
                        round(constraint_upper, 3)
                    )
                    constraint_lower = constraint_upper
                end
                push!(constraints_lower, constraint_lower)
                push!(constraints_upper, constraint_upper)
                push!(present_i, i)
                push!(present_j, j)
            end
        end
    end
    println("Found combined distance constraints from two structures")
    return Constraints(atoms, constraints_lower, constraints_upper, hcat(present_i, present_j))
end


"""
Run the interaction finding pipeline with one or two structures.
If one structure is given, returns a `Constraints` object.
If two structures are given, returns three `Constraints` objects corresponding to the
combined constraints, constraints from structure one and constraints from structure two.
"""
function interactions(pdb_filepath::AbstractString,
                    dssp_filepath::AbstractString;
                    other_ratio::Real=defaults["other_ratio"],
                    tolerance_weight::Real=defaults["tolerance_weight"])
    atoms = readpdb(pdb_filepath)
    dssps = readdssp(dssp_filepath, atoms)
    dists = distances(atoms)
    bs = bonds(atoms, protein_bonds)
    as = angles(atoms, bs)
    divs = divalents(atoms, bs, as)
    inters = findinteractions(atoms, dists, bs, as, divs, dssps)
    return Constraints(atoms, dists, inters, other_ratio=other_ratio, tolerance_weight=tolerance_weight)
end


function interactions(pdb_filepath_one::AbstractString,
                    dssp_filepath_one::AbstractString,
                    pdb_filepath_two::AbstractString,
                    dssp_filepath_two::AbstractString;
                    other_ratio::Real=defaults["other_ratio"],
                    tolerance_weight::Real=defaults["tolerance_weight"])
    atoms_all_one = readpdb(pdb_filepath_one)
    atoms_all_two = readpdb(pdb_filepath_two)
    atoms_one, atoms_two = findcommonatoms(atoms_all_one, atoms_all_two)
    dssps_one = readdssp(dssp_filepath_one, atoms_one)
    dssps_two = readdssp(dssp_filepath_two, atoms_two)
    dists_one = distances(atoms_one)
    bonds_one = bonds(atoms_one, protein_bonds)
    angles_one = angles(atoms_one, bonds_one)
    divalents_one = divalents(atoms_one, bonds_one, angles_one)
    inters_one = findinteractions(atoms_one, dists_one, bonds_one, angles_one, divalents_one, dssps_one)
    constraints_one = Constraints(atoms_one, dists_one, inters_one, other_ratio=other_ratio, tolerance_weight=tolerance_weight)
    dists_two = distances(atoms_two)
    bonds_two = bonds(atoms_two, protein_bonds)
    angles_two = angles(atoms_two, bonds_two)
    divalents_two = divalents(atoms_two, bonds_two, angles_two)
    inters_two = findinteractions(atoms_two, dists_two, bonds_two, angles_two, divalents_two, dssps_two)
    constraints_two = Constraints(atoms_two, dists_two, inters_two, other_ratio=other_ratio, tolerance_weight=tolerance_weight)
    constraints_com = Constraints(atoms_one, dists_one, dists_two, inters_one, inters_two, other_ratio=other_ratio, tolerance_weight=tolerance_weight)
    return constraints_com, constraints_one, constraints_two
end
