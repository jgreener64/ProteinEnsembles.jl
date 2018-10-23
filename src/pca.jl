# Run principal components analysis (PCA)


export
    projectensemble,
    projectstructure,
    projectpdbfile,
    projectpdbfiles,
    fluctuations,
    findimportantpcs


"""
Runs principal components analysis on an ensemble and returns a `PCA` object.
Assumes the ensemble is aligned sensibly.
"""
function PCA(ensemble::ModelledEnsemble)
    # Calculate the covariance matrix
    strucs = ensemble.strucs
    n_strucs = length(strucs)
    @assert n_strucs > 0 "The ensemble does not contain any structures"
    @assert n_strucs > 1 "The ensemble contains 1 structure - need multiple to do PCA"
    av_coords = centroid(ensemble)
    n_atoms = size(av_coords, 2)
    @assert n_atoms == size(strucs[1].coords, 2) "Number of atoms in the ensemble structures and atoms are different"
    inds_to_use = calphaindices(ensemble.atoms)
    n_to_use = length(inds_to_use)
    @assert n_to_use > 0 "Not calculating PCA on any atoms"
    cov = zeros(3*n_to_use, 3*n_to_use)
    av_to_use = getindex(av_coords, collect(1:3), inds_to_use)
    for struc in strucs
        disps = (getindex(struc.coords, collect(1:3), inds_to_use) - av_to_use)[:]
        cov += disps * transpose(disps)
    end
    cov /= n_strucs
    # Find the eigendecomposition
    eig = eigen(cov)
    # Number of non-zero eigenvalues is limited by number of structures
    n_non_zero = min(3*n_to_use-6, n_strucs-1)
    # Reverse the values in order to put the most important PCs first
    pca = PCA(reverse(eig.values)[1:n_non_zero], eig.vectors[:, end:-1:1][:, 1:n_non_zero], av_to_use)
    # Project each structure onto the PCs to get the principal coordinates
    pca.pcs = projectensemble(ensemble, pca)
    return pca
end


"""
Project the structures in an ensemble onto the eigenvectors of a PCA to get the principal coordinates.
Non-zero eigenvalues only are returned.
Assumes the ensemble is aligned and corresponds exactly to the PCA.
"""
function projectensemble(ensemble::ModelledEnsemble, pca::PCA)
    strucs = ensemble.strucs
    n_strucs = length(strucs)
    @assert n_strucs > 0 "The ensemble does not contain any structures"
    av_coords = pca.av_coords
    n_atoms = size(av_coords, 2)
    inds_to_use = calphaindices(ensemble.atoms)
    @assert n_atoms == length(inds_to_use) "Number of atoms in the ensemble and in the PCA are different"
    n_evecs = size(pca.evecs, 2)
    pcs = zeros(n_evecs, n_strucs)
    for j in 1:n_strucs
        disps = getindex(strucs[j].coords, collect(1:3), inds_to_use) - av_coords
        for i in 1:n_evecs
            pcs[i, j] = dot(pca.evecs[:,i], disps[:])
        end
    end
    return pcs
end


"""
Project a list of atoms onto the eigenvectors of a PCA.
Atom map is the mapping of the input atom indices to the reference atoms.
It is not required if the atoms correspond to the PCA.
`calpha_inds` is only required if atom_map is non-trivial and is a list of C-alpha indices.
0 means the atom is not in the reference atoms.
Assumes the structure is aligned.
"""
function projectstructure(atoms::Array{Atom,1},
                    pca::PCA;
                    atom_map::Array{Int,1}=collect(1:length(atoms)),
                    calpha_inds::Array{Int,1}=Int[])
    @assert length(atoms) > 0 "No atoms in atom list"
    # If not all atoms are present, knowledge of the missing c-alphas is required
    if 0 in atom_map
        @assert length(calpha_inds) > 0
    end
    init_to_use = calphaindices(atoms)
    inds_to_use_in = Int[]
    inds_to_use_ref = Int[]
    evec_inds = Int[]
    j = 0
    for i in 1:length(atom_map)
        if atom_map[i] in init_to_use
            j += 1
            push!(inds_to_use_in, atom_map[i])
            push!(inds_to_use_ref, j)
            push!(evec_inds, 3j-2, 3j-1, 3j)
        elseif atom_map[i] == 0 && i in calpha_inds
            j += 1
        end
    end
    n_evecs = size(pca.evecs, 2)
    pcs = zeros(n_evecs, 1)
    coords_to_use = getindex(atomcoords(atoms), collect(1:3), inds_to_use_in)
    disps = coords_to_use - getindex(pca.av_coords, collect(1:3), inds_to_use_ref)
    for i in 1:n_evecs
        pcs[i] = dot(getindex(pca.evecs, evec_inds, i), disps[:])
    end
    return pcs
end


"""
Projects the structure in a PDB file onto the eigenvectors of a PCA.
The residue numbering between the PDB files must be consistent but missing/extra
residues are allowed.
Reads the PDB file, maps the atoms to the ref atoms, aligns the atoms to the ref
atoms then gets the projections.
Arguments are PDB filepath, reference atoms to get atom map, coords to align to
and the PCA object.
"""
function projectpdbfile(pdb_filepath::AbstractString,
                    atoms_ref::Array{Atom,1},
                    ens_al::Array{Float64},
                    pca::PCA)
    @assert isfile(pdb_filepath) "Not a valid filepath: \"$pdb_filepath\""
    atoms_in = readpdb(pdb_filepath)
    atom_map = atommap(atoms_in, atoms_ref)
    # Input atoms not in the ref atoms are removed to allow alignment
    atoms_in_min = Atom[]
    ens_al_min = Float64[]
    atom_map_new = Int[]
    j = 0
    for i in 1:length(atoms_ref)
        if atom_map[i] != 0
            push!(atoms_in_min, atoms_in[atom_map[i]])
            push!(ens_al_min, ens_al[:,i]...)
            j += 1
            push!(atom_map_new, j)
        else
            push!(atom_map_new, 0)
        end
    end
    alignatoms!(atoms_in_min, reshape(ens_al_min, 3, j))
    pcs = projectstructure(atoms_in_min, pca, atom_map=atom_map_new, calpha_inds=calphaindices(atoms_ref))
    return pcs
end


"Wrapper function for projectpdbfile to project multiple PDBs and write the output."
function projectpdbfiles(out_prefix::AbstractString,
                    pdb_filepaths::Array{<:AbstractString,1},
                    atoms_ref::Array{Atom,1},
                    ens_al::Array{Float64},
                    pca::PCA)
    checkfilepath(out_prefix)
    @assert length(pdb_filepaths) > 0 "No entries in pdb_filepaths"
    pcs_all = Float64[]
    for (i, pdb_filepath) in enumerate(pdb_filepaths)
        pcs = projectpdbfile(pdb_filepath, atoms_ref, ens_al, pca)
        writeprojections("$(out_prefix)$i.tsv", pcs)
        if length(pcs_all) > 0
            pcs_all = hcat(pcs_all, pcs)
        else
            pcs_all = pcs[:]
        end
    end
    return pcs_all
end


"""
Gets mean square fluctuation for each atom in an ensemble.
Assumes the ensemble is aligned sensibly.
Returns a list of mean square fluctuations.
"""
function fluctuations(ensemble::ModelledEnsemble)
    strucs = ensemble.strucs
    n_strucs = length(strucs)
    @assert n_strucs > 0 "The ensemble does not contain any structures"
    @assert n_strucs > 1 "The ensemble contains 1 structure - need multiple to get fluctuations"
    av_coords = centroid(ensemble)
    n_atoms = size(av_coords, 2)
    @assert n_atoms == size(strucs[1].coords, 2) "Number of atoms in the ensemble and in the PCA are different"
    inds_to_use = calphaindices(ensemble.atoms)
    n_to_use = length(inds_to_use)
    @assert n_to_use > 0 "Not calculating fluctuations on any atoms"
    flucs = zeros(1, n_to_use)
    av_to_use = getindex(av_coords, collect(1:3), inds_to_use)
    for struc in strucs
        disps = (getindex(struc.coords, collect(1:3), inds_to_use) - av_to_use)
        flucs += sum(disps .* disps, dims=1)
    end
    flucs /= n_strucs
    return flucs[:]
end


"""
Takes two sets of projections onto the PCs of a PCA and returns the PC indices ordered
such that the first represents the largest difference between the two structures.
"""
function findimportantpcs(pcs_one::Array{Float64}, pcs_two::Array{Float64})
    @assert size(pcs_one) == size(pcs_two) "PC arrays have different sizes"
    @assert size(pcs_one, 2) == 1 "Each set of projections should be of one structure"
    diff = abs.(pcs_one - pcs_two)
    pc_ord = sortperm(diff[:], rev=true)
    return pc_ord
end
