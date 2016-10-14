# Align structures and ensembles


export
    kabschalignment,
    displacements,
    align!,
    alignsimple!,
    alignatoms!,
    alignensemble!,
    alignensemblesimple!,
    centroid,
    rmsd,
    ensemblermsds,
    selfalignensemble!


"""
Find the transformations to map the first set of coordinates onto the second.
See Kabsch algorithm for details.
Returns the translation to center the first set, the translation to center the second set
and the rotation to map the first set onto the second.
"""
function kabschalignment(coords_one::Array{Float64}, coords_two::Array{Float64})
    n_coords = size(coords_one, 2)
    @assert n_coords == size(coords_two, 2) "Number of coordinates is different"
    @assert n_coords > 0 "The coordinate arrays are empty"
    # Find the translation to move the coordinate centroids to the origin
    trans_one = sum(coords_one, 2) / n_coords
    trans_two = sum(coords_two, 2) / n_coords
    p = coords_one - repeat(trans_one, inner=[1,n_coords])
    q = coords_two - repeat(trans_two, inner=[1,n_coords])
    # Find the rotation that maps the coordinates
    cov = p * transpose(q)
    svd = svdfact(cov)
    rotation = svd[:V] * transpose(svd[:U])
    return trans_one, trans_two, rotation
end


"Return the distances between each point for two sets of aligned coordinates."
function displacements(coords_one::Array{Float64}, coords_two::Array{Float64})
    n_coords = size(coords_one, 2)
    @assert n_coords == size(coords_two, 2) "Number of coordinates is different"
    @assert n_coords > 0 "The coordinate arrays are empty"
    diff = coords_one - coords_two
    sq_diff = diff .* diff
    dists = sqrt(sum(sq_diff, 1))
    return dists[:]
end


"""
Align the first set of atomic coordinates to the second using an iterative procedure.
At each step atoms with a high structural deviation are removed.
If the standard deviation of the structural deviations is greater than `cutoff`, then structural
deviations further than 2 standard deviations from the median are discarded.
This is repeated `n_cycles` times.
If `calpha` is true, only considers the C-alpha atoms.
"""
function align!(coords_one::Array{Float64},
                    coords_two::Array{Float64},
                    atoms::Array{Atom,1};
                    calpha::Bool=true,
                    n_cycles::Integer=defaults["align_cycles"],
                    cutoff::Real=defaults["align_cutoff"])
    n_coords = size(coords_one, 2)
    @assert n_coords == size(coords_two, 2) "Number of coordinates is different"
    @assert n_coords > 0 "The coordinate arrays are empty"
    n_atoms = length(atoms)
    @assert n_coords == n_atoms "Number of coordinates and number of atoms are different"
    if calpha
        inds_to_use = calphaindices(atoms)
    else
        inds_to_use = collect(1:n_atoms)
    end
    n_starting = length(inds_to_use)
    @assert n_starting > 0 "Not aligning on any atoms"
    new_coords = zeros(3, n_atoms)
    for i in 1:n_cycles
        # Find the alignment
        trans_one, trans_two, rot = kabschalignment(getindex(coords_one, collect(1:3), inds_to_use), getindex(coords_two, collect(1:3), inds_to_use))
        # Move the first set to the origin, rotate it then move it to the centroid of the second set
        new_coords = coords_one - repeat(trans_one, inner=[1,n_atoms])
        new_coords = rot * new_coords + repeat(trans_two, inner=[1,n_atoms])
        #println("Iteration ", i, " RMSD: ", round(rmsd(getindex(new_coords, collect(1:3), inds_to_use), getindex(coords_two, collect(1:3), inds_to_use)), 1), " A")
        struc_devs = displacements(getindex(new_coords, collect(1:3), inds_to_use), getindex(coords_two, collect(1:3), inds_to_use))
        sd = sqrt(var(struc_devs))
        dev_cutoff = median(struc_devs) + 2*sd
        inds_to_remove = Int[]
        for j in 1:length(struc_devs)
            if struc_devs[j] > dev_cutoff
                push!(inds_to_remove, j)
            end
        end
        if length(inds_to_remove) > 0 && sd > cutoff && i != n_cycles
            #println("Removing ", length(inds_to_remove), " atoms")
            deleteat!(inds_to_use, inds_to_remove)
            @assert length(inds_to_use) > 0 "Not aligning on any atoms"
        elseif i == n_cycles
            #println("Finishing after ", n_cycles, " iterations")
        else
            #println("Stopping early")
            break
        end
    end
    #println("Alignment has RMSD ", round(rmsd(getindex(new_coords, collect(1:3), inds_to_use), getindex(coords_two, collect(1:3), inds_to_use)), 1), " A based on ", length(inds_to_use), " out of ", n_starting, " atoms")
    coords_one[:] = new_coords[:]
end


"""
Align the first set of atomic coordinates to the second using a single Kabsch rotation.
If `calpha` is true, only considers the C-alpha atoms.
"""
function alignsimple!(coords_one::Array{Float64},
                    coords_two::Array{Float64},
                    atoms::Array{Atom,1};
                    calpha::Bool=true)
    n_coords = size(coords_one, 2)
    @assert n_coords == size(coords_two, 2) "Number of coordinates is different"
    @assert n_coords > 0 "The coordinate arrays are empty"
    n_atoms = length(atoms)
    @assert n_coords == n_atoms "Number of coordinates and number of atoms are different"
    if calpha
        inds_to_use = calphaindices(atoms)
    else
        inds_to_use = collect(1:n_atoms)
    end
    @assert length(inds_to_use) > 0 "Not aligning on any atoms"
    # Find the alignment
    trans_one, trans_two, rot = kabschalignment(getindex(coords_one, collect(1:3), inds_to_use), getindex(coords_two, collect(1:3), inds_to_use))
    # Move the first set to the origin, rotate it then move it to the centroid of the second set
    new_coords = coords_one - repeat(trans_one, inner=[1,n_atoms])
    new_coords = rot * new_coords + repeat(trans_two, inner=[1,n_atoms])
    coords_one[:] = new_coords[:]
end


"""
Align an array of atoms to a set of coordinates or another set of atoms.
Input is list of atoms and coords to align to.
Changes the atom coords.
"""
function alignatoms!(atoms::Array{Atom,1}, coords_ref::Array{Float64})
    n_atoms = length(atoms)
    @assert n_atoms > 0 "There are no atoms in the atom list"
    @assert n_atoms == size(coords_ref, 2) "The number of atoms and reference coordinates differ"
    coords = atomcoords(atoms)
    align!(coords, coords_ref, atoms)
    atomcoords!(atoms, coords)
end


alignatoms!(atoms::Array{Atom,1}, atoms_ref::Array{Atom,1}) = alignatoms!(atoms, atomcoords(atoms_ref))


"Align an ensemble of coordinates to a coordinate set using an iterative procedure for each structure."
function alignensemble!(ensemble::ModelledEnsemble, coords_ref::Array{Float64})
    strucs = ensemble.strucs
    @assert length(strucs) > 0 "No structures in the ensemble"
    n_coords_ref = size(coords_ref, 2)
    @assert n_coords_ref > 0 "No coords in reference coordinate array"
    @assert size(strucs[1].coords, 2) == n_coords_ref "Number of atoms in ensemble and reference structure are different"
    atoms = ensemble.atoms
    for struc in strucs
        align!(struc.coords, coords_ref, atoms)
    end
end


"Align an ensemble of coordinates to a coordinate set using a single Kabsch rotation for each structure."
function alignensemblesimple!(ensemble::ModelledEnsemble, coords_ref::Array{Float64})
    strucs = ensemble.strucs
    @assert length(strucs) > 0 "No structures in the ensemble"
    n_coords_ref = size(coords_ref, 2)
    @assert n_coords_ref > 0 "No coords in reference coordinate array"
    @assert size(strucs[1].coords, 2) == n_coords_ref "Number of atoms different in ensemble and reference structure"
    atoms = ensemble.atoms
    for struc in strucs
        alignsimple!(struc.coords, coords_ref, atoms)
    end
end


"""
Find the average structure in an ensemble.
The ensemble should be aligned.
"""
function centroid(ensemble::ModelledEnsemble)
    strucs = ensemble.strucs
    @assert length(strucs) > 0 "No structures in the ensemble"
    n_coords = size(strucs[1].coords, 2)
    @assert n_coords > 0 "No coordinates for one or more of the structures"
    sum_coords = zeros(3, n_coords)
    for struc in strucs
        @assert n_coords == size(struc.coords, 2) "There are varying numbers of coordinates between structures"
        sum_coords += struc.coords
    end
    sum_coords /= length(strucs)
    return sum_coords
end


"""
Get the RMSD between two identical structures or two sets of coordinates.
Assumes the structures are aligned.
If `calpha` is true, only considers the C-alpha atoms.
"""
function rmsd(coords_one::Array{Float64},
                    coords_two::Array{Float64},
                    atoms::Array{Atom,1};
                    calpha::Bool=true)
    n_coords = size(coords_one, 2)
    @assert n_coords == size(coords_two, 2) "Number of coordinates is different"
    @assert n_coords > 0 "The coordinate arrays are empty"
    n_atoms = length(atoms)
    @assert n_coords == n_atoms "Number of coordinates and number of atoms are different"
    if calpha
        inds_to_use = calphaindices(atoms)
    else
        inds_to_use = collect(1:n_atoms)
    end
    return rmsd(getindex(coords_one, collect(1:3), inds_to_use), getindex(coords_two, collect(1:3), inds_to_use))
end

function rmsd(coords_one::Array{Float64}, coords_two::Array{Float64})
    n_coords = size(coords_one, 2)
    @assert n_coords == size(coords_two, 2) "Number of coordinates is different"
    @assert n_coords > 0 "The coordinate arrays are empty"
    diff = coords_one - coords_two
    sq_diff = diff .* diff
    return sqrt(sum(sq_diff) / n_coords)
end


"""
Get the RMSD between each structure in an ensemble and a reference structure.
"""
function ensemblermsds(ensemble::ModelledEnsemble, atoms_ref::Array{Atom,1})
    @assert length(atoms_ref) == length(ensemble.atoms) "Number of atoms in ensemble and reference are different"
    rmsds = Float64[]
    coords_ref = atomcoords(atoms_ref)
    for struc in ensemble.strucs
        push!(rmsds, rmsd(struc.coords, coords_ref, atoms_ref))
    end
    return rmsds
end


"""
Align an ensemble to itself.
See Bakan and Bahar, PNAS, 2009 for algorithm.
The single Kabsch alignment is used during iteration, then the iterative alignment is done at the end.
Returns the structure the final alignment was on, not the ensemble average (but they will be close).
The cutoff parameter determines the RMSD between two average structures below which iteration is terminated.
"""
function selfalignensemble!(ensemble::ModelledEnsemble; cutoff::Real=defaults["ens_align_cutoff"])
    strucs = ensemble.strucs
    n_strucs = length(strucs)
    @assert n_strucs > 0 "The ensemble does not contain any structures"
    @assert n_strucs > 1 "The ensemble contains 1 structure - need multiple for alignment"
    atoms = ensemble.atoms
    @assert length(atoms) > 0 "The atom list for the ensemble is empty"
    chosen_struc = strucs[rand(1:n_strucs)]
    coords_old = chosen_struc.coords
    alignensemblesimple!(ensemble, coords_old)
    current_rmsd = cutoff+1.0 # Initialised with arbitrary value greater than the cutoff
    counter = 0
    coords_new = Float64[]
    while current_rmsd >= cutoff
        if counter > 1000
            error("Could not align ensemble to itself")
        end
        coords_new = centroid(ensemble)
        alignensemblesimple!(ensemble, coords_new)
        current_rmsd = rmsd(coords_new, coords_old, atoms)
        coords_old = coords_new
        counter += 1
    end
    alignensemble!(ensemble, coords_new)
    println("Aligned ensemble to itself after ", counter, " iteration(s)")
    return coords_new
end
