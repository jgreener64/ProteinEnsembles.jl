# Generate new interactions to perturb ensembles


export
    repeatpocketpoints,
    perturbensemble


"""
Find new interactions arising from fake modulator atoms and form a new `Constraints` object.
Arguments are atoms to calculate constraints with, `Constraints` object to add to, generated modulator
coordinates, distance cutoff and constraint tolerances for intra-modulator and modulator-protein
interactions, and minimum lower distance constraint.
Returns a new `Constraints` object including the modulator.
"""
function Constraints(atoms::Array{Atom,1},
                    old_constraints::Constraints,
                    mod_coords::Array{Float64};
                    intra_cutoff::Real=defaults["mod_intra_cutoff"], # If this is 0.0 there are no intra-modulator interactions
                    intra_tolerance::Real=defaults["mod_intra_tolerance"],
                    inter_cutoff::Real=defaults["mod_inter_cutoff"],
                    inter_tolerance::Real=defaults["mod_inter_tolerance"],
                    min_constraint_dist::Real=defaults["mod_min_constraint_dist"])
    n_atoms = length(atoms)
    @assert n_atoms > 0 "No atoms in atom list"
    constraint_atoms = old_constraints.atoms
    @assert n_atoms == length(constraint_atoms) "Different number of atoms in atom list and constraints"
    n_constraints = size(old_constraints.pres_inds, 1)
    @assert n_constraints == length(old_constraints.lower) && n_constraints == length(old_constraints.upper) "The number of lower and upper constraints do not correspond to the number of constraints present"
    n_mod_atoms = size(mod_coords, 2)
    @assert n_mod_atoms > 0 "No modulator coordinates in list"
    @assert intra_cutoff >= 0 "intra_cutoff cannot be negative"
    @assert intra_tolerance >= 0 "intra_tolerance cannot be negative"
    @assert inter_cutoff >= 0 "inter_cutoff cannot be negative"
    @assert inter_tolerance >= 0 "inter_tolerance cannot be negative"
    @assert min_constraint_dist >= 0 "min_constraint_dist cannot be negative"
    n_total = n_atoms + n_mod_atoms
    # Copy Constraints object and form larger constraints lists
    new_lower = copy(old_constraints.lower)
    new_upper = copy(old_constraints.upper)
    new_pres_inds = copy(old_constraints.pres_inds)
    new_atoms = deepcopy(constraint_atoms)
    # Add new atoms to atom array
    for i in 1:n_mod_atoms
        push!(new_atoms, Atom(mod_atom_info["atom_name"], mod_atom_info["res_name"], mod_atom_info["chain_id"], mod_atom_info["res_n"], mod_coords[:,i], mod_atom_info["element"]))
    end
    new_inds_i = Int[]
    new_inds_j = Int[]
    # Find interactions within modulator
    intra_counter = 0
    intra_cutoff_sq = intra_cutoff^2
    for i in 1:n_mod_atoms
        for j in 1:(i-1)
            sq_dist = 0.0
            for m = 1:3
                sq_dist += abs2(mod_coords[m,i] - mod_coords[m,j])
            end
            if sq_dist < intra_cutoff_sq
                dist = sqrt(sq_dist)
                push!(new_lower, max(dist-intra_tolerance, min_constraint_dist))
                push!(new_upper, dist+intra_tolerance)
                push!(new_inds_i, n_atoms+i)
                push!(new_inds_j, n_atoms+j)
                intra_counter += 1
            end
        end
    end
    # Find interactions between modulator and protein
    inter_counter = 0
    inter_cutoff_sq = inter_cutoff^2
    for i in 1:n_mod_atoms
        for j in 1:n_atoms
            sq_dist = 0.0
            for m = 1:3
                sq_dist += abs2(mod_coords[m,i] - atoms[j].coords[m])
            end
            if sq_dist < inter_cutoff_sq
                dist = sqrt(sq_dist)
                push!(new_lower, max(dist-inter_tolerance, min_constraint_dist))
                push!(new_upper, dist+inter_tolerance)
                push!(new_inds_i, n_atoms+i)
                push!(new_inds_j, j)
                inter_counter += 1
            end
        end
    end
    if length(new_pres_inds) > 0
        new_pres_inds = vcat(new_pres_inds, hcat(new_inds_i, new_inds_j))
    else
        new_pres_inds = hcat(new_inds_i, new_inds_j)
    end
    #println("Found ", intra_counter, " interactions within the modulator")
    println("Found ", inter_counter, " interactions between the modulator and the protein")
    return Constraints(new_atoms, new_lower, new_upper, new_pres_inds)
end


"""
Repeat pocket point to get a certain number of points.
Points are repeated as long as whole copies can be made, then the remainder are chosen randomly.
Returns coordinate array.
"""
function repeatpocketpoints(coords::Array{Float64}, n_out_points::Integer=defaults["mod_n_points"])
    n_pocket_points = Int(length(coords) / 3)
    inds_to_use = repeat(collect(1:n_pocket_points), outer=[Int(floor(n_out_points / n_pocket_points))])
    append!(inds_to_use, rand(1:n_pocket_points, n_out_points % n_pocket_points))
    return getindex(coords, collect(1:3), inds_to_use)
end


"""
Gets ensembles for multiple modulators.
Use pocket points to generate additional constraints.
"""
function perturbensemble(atoms::Array{Atom,1},
                    constraints::Constraints,
                    n_strucs::Integer,
                    mod_path::Union{AbstractString, Void},
                    n_mods::Integer)
    ensemble_mods = ModelledEnsemble[]
    if n_mods > 0 && mod_path != nothing
        pock_points = readpocketpoints(mod_path)
        n_pocks = length(pock_points)
        if n_mods > n_pocks
            println(n_mods, " modulators asked for but only ", n_pocks, " found")
        end
        n_mods_to_use = min(n_mods, n_pocks)
        for i in 1:n_mods_to_use
            println("Generating ensemble with modulator ", i, " of ", n_mods_to_use)
            mod_coords = repeatpocketpoints(pock_points[i])
            new_constraints = Constraints(atoms, constraints, mod_coords)
            push!(ensemble_mods, generateensemble(new_constraints, n_strucs))
        end
    end
    return ensemble_mods
end
