# All types used


export
    BondedPair,
    Atom,
    Constraints,
    ModelledStructure,
    ModelledEnsemble,
    PCA


"A pair of atoms with a covalent bond between them in the standard amino acids."
type BondedPair
    residue::String
    atom_one::String
    atom_two::String
end

function Base.show(io::IO, bonded_pair::BondedPair)
    println(io, "Bond between $(bonded_pair.atom_one) and $(bonded_pair.atom_two) in $(bonded_pair.residue)")
end


"An atom from a PDB format file."
type Atom
    atom_name::String
    res_name::String
    chain_id::Char
    res_n::Int
    coords::Array{Float64,1}
    element::String
end

function Base.show(io::IO, atom::Atom)
    println(io, "Atom name      - ", atom.atom_name)
    println(io, "Residue name   - ", atom.res_name)
    println(io, "Chain ID       - ", atom.chain_id)
    println(io, "Residue number - ", atom.res_n)
    println(io, "Coordinates    - [", atom.coords[1], ", ", atom.coords[2], ", ", atom.coords[3], "]")
    println(io, "Element        - ", atom.element)
end


"Information about the lower and upper distance constraints between atoms."
type Constraints
    atoms::Array{Atom,1}
    lower::Array{Float64,1} # List of lower constraints by present index
    upper::Array{Float64,1} # List of upper constraints by present index
    pres_inds::Array{Int,2} # List of index pairs of present constraints
end

function Base.show(io::IO, constraints::Constraints)
    println(io, length(constraints.lower), " lower and upper distance constraints for ", length(constraints.atoms), " atoms")
end


"A set of coordinates and a score."
type ModelledStructure
    score::Float64 # This is sometimes set as a placeholder of -1.0 to mean no score
    coords::Array{Float64}
end

function Base.show(io::IO, struc::ModelledStructure)
    println(io, "Modelled structure with score ", round(struc.score, 1), " and a ", size(struc.coords), " coordinate array")
end


"A group of modelled structures."
type ModelledEnsemble
    atoms::Array{Atom,1}
    strucs::Array{ModelledStructure,1}
end

# Initialise with an empty list of structures
ModelledEnsemble(atoms::Array{Atom,1}) = ModelledEnsemble(atoms, ModelledStructure[])

function Base.show(io::IO, ens::ModelledEnsemble)
    println(io, "Ensemble of ", length(ens.strucs), " structures from ", length(ens.atoms), " atoms")
end


"The results of a principal components analysis, sorted with highest eigenvalue first."
type PCA
    evals::Array{Float64}
    evecs::Array{Float64}
    av_coords::Array{Float64} # Average C-alpha coordinates
    pcs::Array{Float64} # The principal coordinates
end

# Initialise with an empty list of PCs
PCA(evals::Array{Float64}, evecs::Array{Float64}, av_coords::Array{Float64}) = PCA(
    evals, evecs, av_coords, Float64[])

function Base.show(io::IO, pca::PCA)
    println(io, "PCA with ", length(pca.evals), " principal components for ", size(pca.av_coords, 2), " atoms across ", size(pca.pcs, 2), " structures")
end
