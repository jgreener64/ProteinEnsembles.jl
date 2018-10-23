# Generate structures that fill distance constraints


export
    fillconstraints,
    errorscore,
    fixchirality!,
    addstructure!,
    removemodulator!,
    generateensemble,
    trimbyscore!


using ProgressMeter


"""
Use stochastic proximity embedding to generate a structure iteratively.
See Agrafiotis et.al., Chapter 14 in Mucherino et. al., 2013 for algorithm details (Algorithm 8).
Returns coordinate list.
"""
function fillconstraints(constraints::Constraints;
                    iters_per_atom::Integer=defaults["iters_per_atom"],
                    cyc_step_ratio::Real=defaults["cyc_step_ratio"],
                    box_size::Real=defaults["box_size"]
    )
    @assert iters_per_atom > 0 "iters_per_atom must be positive"
    @assert cyc_step_ratio > 0 "cyc_step_ratio must be positive"
    constraints_lower = constraints.lower
    constraints_upper = constraints.upper
    present_i = constraints.pres_inds[:,1]
    present_j = constraints.pres_inds[:,2]
    squares_lower = constraints_lower .* constraints_lower
    squares_upper = constraints_upper .* constraints_upper
    n_atoms = length(constraints.atoms)
    # Work out n_cycles and n_steps using n_atoms, iters_per_atom and the ratio of steps to cycles
    n_cycs_steps = n_atoms * iters_per_atom
    n_cycles = round(Int, sqrt(n_cycs_steps / cyc_step_ratio))
    n_steps = round(Int, sqrt(n_cycs_steps * cyc_step_ratio))
    # Randomise atomic coordinates in a cube
    coords = rand!(zeros(3, n_atoms)) * box_size
    n_constraints = length(present_i)
    learning_rate = 1.0
    coord_new_i = zeros(3)
    coord_new_j = zeros(3)
    constraint_near = 0.0
    rand_list = zeros(Int, n_steps)
    @inbounds for cyc in 1:n_cycles
        # List of random index tuples
        rand!(rand_list, 1:n_constraints)
        for rand_ind in rand_list
            # Pick a random constraint
            i = present_i[rand_ind]
            j = present_j[rand_ind]
            # Compute distance between atoms
            sq_dist = 0.0
            for m = 1:3
                sq_dist += abs2(coords[m, i] - coords[m, j])
            end
            # Check if distance is within constraints
            if sq_dist < squares_lower[rand_ind]
                constraint_near = constraints_lower[rand_ind]
            elseif sq_dist > squares_upper[rand_ind]
                constraint_near = constraints_upper[rand_ind]
            else
                continue
            end
            # Adjust coordinates
            dist = sqrt(sq_dist)
            prefactor = learning_rate * 0.5 * (constraint_near - dist) / dist
            for k in 1:3
                coord_new_i[k] = newcoords(k, i, j, prefactor, coords)
                coord_new_j[k] = newcoords(k, j, i, prefactor, coords)
            end
            for l in 1:3
                coords[l, i] = coord_new_i[l]
                coords[l, j] = coord_new_j[l]
            end
        end
        # Decrease learning rate
        learning_rate -= 1 / n_cycles
    end
    return coords
end


# Speeds up fillconstraints
newcoords(k::Integer, i::Integer, j::Integer, prefactor::Float64, coords::Array{Float64,2}) = coords[k, i] + prefactor * (coords[k, i] - coords[k, j])


"Calculate SPE error score."
function errorscore(constraints::Constraints, coords::Array{Float64})
    constraints_lower = constraints.lower
    constraints_upper = constraints.upper
    pres_inds = constraints.pres_inds
    score = 0.0
    for m in 1:size(pres_inds, 1)
        k, l = pres_inds[m,:]
        dist = norm(coords[:,k] - coords[:,l])
        if dist < constraints_lower[m]
            # Arbitrary minimum constraint size used here to prevent division by zero and infinite scores
            score += ((dist-constraints_lower[m])^2) / max((constraints_upper[m]-constraints_lower[m]), 0.001)
        elseif dist > constraints_upper[m]
            score += ((dist-constraints_upper[m])^2) / max((constraints_upper[m]-constraints_lower[m]), 0.001)
        end
    end
    return score
end


"""
Checks the chirality of coordinate set and switches it if it isn't correct.
The coordinates must correspond to the atoms.
Inverts each z coordinate in the coordinate array if required.
"""
function fixchirality!(coords::Array{Float64}, atoms::Array{Atom,1})
    n_atoms = length(atoms)
    @assert n_atoms > 0 "No atoms in atom list"
    @assert n_atoms == size(coords, 2) "Number of atoms and coordinates are different"
    count_l = 0
    count_d = 0
    for outer_id in 1:n_atoms
        if atoms[outer_id].atom_name == "CA"
            res_n = atoms[outer_id].res_n
            chain_id = atoms[outer_id].chain_id
            coord_ca = coords[:, outer_id]
            found_co = false
            found_n = false
            found_cb = false
            vec_oc = zeros(3)
            vec_nc = zeros(3)
            vec_rc = zeros(3)
            for inner_id in 1:n_atoms
                if atoms[inner_id].res_n == res_n && atoms[inner_id].chain_id == chain_id
                    if atoms[inner_id].atom_name == "C"
                        found_co = true
                        vec_oc = coord_ca - coords[:,inner_id]
                    elseif atoms[inner_id].atom_name == "N"
                        found_n = true
                        vec_nc = coord_ca - coords[:,inner_id]
                    elseif atoms[inner_id].atom_name == "CB"
                        found_cb = true
                        vec_rc = coord_ca - coords[:,inner_id]
                    end
                end
            end
            if found_co && found_n && found_cb
                if dot(vec_oc, LinearAlgebra.cross(vec_rc, vec_nc)) > 0
                    count_d += 1
                else
                    count_l += 1
                end
            end
        end
    end
    @assert max(count_d, count_l) > 0 "Could not find any tetrahedral C-alpha systems to determine chirality"
    if count_d > count_l
        coords[3,:] *= -1
    end
end


"Add a structure to an ensemble."
function addstructure!(ensemble::ModelledEnsemble, struc::ModelledStructure)
    n_coords = size(struc.coords, 2)
    @assert n_coords > 0 "The structure to add does not contain any coordinates"
    strucs = ensemble.strucs
    if length(strucs) > 0
        @assert n_coords == size(strucs[1].coords, 2) "Attempting to add a structure with a different number of atoms to the existing structures"
    end
    push!(strucs, struc)
end


"Takes a `ModelledEnsemble` and removes the fake modulator atoms and coords."
function removemodulator!(ensemble::ModelledEnsemble)
    atoms = ensemble.atoms
    n_atoms = length(atoms)
    @assert n_atoms > 0 "There are no atoms in the atom list"
    strucs = ensemble.strucs
    n_strucs = length(strucs)
    @assert n_strucs == 0 || size(strucs[1].coords, 2) == n_atoms "Number of coordinates and number of atoms are different"
    # Find indices corresponding to non-modulator atoms
    atom_inds = Int[]
    for i in 1:n_atoms
        if !ismodulator(atoms[i])
            push!(atom_inds, i)
        end
    end
    if length(atom_inds) != n_atoms
        # Remove modulator atoms from atom list
        ensemble.atoms = getindex(atoms, atom_inds)
        # Remove modulator coords from each structure
        for j in 1:n_strucs
            ensemble.strucs[j].coords = getindex(strucs[j].coords, collect(1:3), atom_inds)
        end
    end
end


"""
Generate sets of coordinates with the correct chirality from distance constraints.
discard_ratio gives the fraction of extra structures to calculate, then the low-scoring structures are discarded.
If `remove_mod` is true then modulator atoms and coords are removed from the ensemble before it is returned.
Returns a `ModelledEnsemble`.
"""
function generateensemble(constraints::Constraints,
                    n_strucs::Integer;
                    iters_per_atom::Integer=defaults["iters_per_atom"],
                    discard_ratio::Real=defaults["discard_ratio"],
                    remove_mod::Bool=true)
    atoms = constraints.atoms
    n_atoms = length(atoms)
    @assert n_atoms > 0 "The Constraints object does not have the atoms set"
    n_constraints = size(constraints.pres_inds, 1)
    @assert n_constraints == length(constraints.lower) && n_constraints == length(constraints.upper) "The number of lower and upper constraints do not correspond to the number of constraints present"
    @assert n_constraints > 0 "There are no constraints present"
    @assert discard_ratio >= 1.0 "discard_ratio must be >= 1"
    ensemble = ModelledEnsemble(atoms)
    n_strucs_init = round(Int, ceil(n_strucs * discard_ratio))
    println("Will generate ", n_strucs_init, " structures and keep best-scoring ", n_strucs)
    p = Progress(n_strucs_init, 1, "Progress: ", 50)
    for struc_n in 1:n_strucs_init
        coords = fillconstraints(constraints, iters_per_atom=iters_per_atom)
        score = errorscore(constraints, coords)
        fixchirality!(coords, atoms)
        addstructure!(ensemble, ModelledStructure(score, coords))
        next!(p)
    end
    trimbyscore!(ensemble, n_strucs)
    if remove_mod
        removemodulator!(ensemble)
    end
    return ensemble
end


"Retains only the n_to_keep top scoring coordinate sets in an ensemble."
function trimbyscore!(ensemble::ModelledEnsemble, n_to_keep::Integer)
    strucs = ensemble.strucs
    @assert n_to_keep > 0 "n_to_keep must be positive"
    @assert n_to_keep <= length(strucs) "Asking to trim an ensemble of size $(length(strucs)) to $n_to_keep structures"
    scores = [struc.score for struc in strucs]
    inds_to_keep = sortperm(scores)[1:n_to_keep]
    ensemble.strucs = getindex(strucs, inds_to_keep)
    println("Kept lowest scoring ", n_to_keep, " out of ", length(scores), " structures")
    println("Scores range from ", round(sort(scores)[1], digits=1), " to ", round(sort(scores)[n_to_keep], digits=1))
end
