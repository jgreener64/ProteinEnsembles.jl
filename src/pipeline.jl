# Pipeline wrapper functions


export
    runfromshell,
    paramfromshell,
    runpipeline,
    runanalysis


using ArgParse


const exprose_repo_url = "https://github.com/jgreener64/ProteinEnsembles.jl"


"Read arguments from command line for an ExProSE run."
function parsecommandlinerun()
    s = ArgParseSettings(prog="julia run.jl",
                    description="Run the ExProSE procedure to generate ensembles of protein structures. See $exprose_repo_url for the ExProSE docs and citation.")
    @add_arg_table s begin
        "--i1"
            help = "filepath to the first input PDB file (required)"
            arg_type = AbstractString
            required = true
        "--d1"
            help = "filepath to the DSSP file for the first input PDB file (required)"
            arg_type = AbstractString
            required = true
        "--i2"
            help = "filepath to the second input PDB file (optional)"
            arg_type = AbstractString
        "--d2"
            help = "filepath to the DSSP file for the second input PDB file (optional)"
            arg_type = AbstractString
        "--out_dir", "-o"
            help = "directory to write output to"
            arg_type = AbstractString
            default = defaults["out_dir"]
        "--n_strucs", "-n"
            help = "number of structures to generate"
            arg_type = Int
            default = defaults["n_strucs"]
        "--tolerance_weight", "-w"
            help = "weighting of constraint tolerances for interactions"
            arg_type = Float64
            default = defaults["tolerance_weight"]
        "--other_ratio", "-r"
            help = "ratio of non-specific interactions to atom number"
            arg_type = Float64
            default = defaults["other_ratio"]
        "--extra_pdbs", "-e"
            help = "additional PDB files to be projected onto the PCs of the ensemble"
            arg_type = AbstractString
            nargs = '+'
        "--mod_path", "-l"
            help = "path to PDB file for modulator point locations; set this, along with n_mods, to perturb the ensemble"
            arg_type = AbstractString
        "--n_mods", "-m"
            help = "number of modulators to separately perturb the ensemble with; set this, along with mod_path, to perturb the ensemble"
            arg_type = Int
            default = 0
    end
    return parse_args(s)
end


"Read arguments from command line for an ExProSE auto-parameterisation."
function parsecommandlineparam()
    s = ArgParseSettings(prog="julia param.jl",
                    description="Run the ExProSE auto-parameterisation procedure. See $exprose_repo_url for the ExProSE docs and citation.")
    @add_arg_table s begin
        "--i1"
            help = "filepath to the first input PDB file (required)"
            arg_type = AbstractString
            required = true
        "--d1"
            help = "filepath to the DSSP file for the first input PDB file (required)"
            arg_type = AbstractString
            required = true
        "--i2"
            help = "filepath to the second input PDB file (required)"
            arg_type = AbstractString
            required = true
        "--d2"
            help = "filepath to the DSSP file for the second input PDB file (required)"
            arg_type = AbstractString
            required = true
        "--out_dir", "-o"
            help = "directory to write output to"
            arg_type = AbstractString
            default = defaults["out_dir_param"]
        "--n_strucs", "-n"
            help = "number of structures to generate"
            arg_type = Int
            default = defaults["n_strucs_param"]
        "--other_ratio", "-r"
            help = "ratio of non-specific interactions to atom number"
            arg_type = Float64
            default = defaults["other_ratio"]
        "--fraction", "-f"
            help = "fraction of structures between input structures to search for during parameterisation"
            arg_type = Float64
            default = defaults["frac_between"]
        "--tmscore", "-t"
            help = "executable to run TMscore, see http://zhanglab.ccmb.med.umich.edu/TM-score"
            arg_type = AbstractString
            default = defaults["tmscore_path"]
    end
    return parse_args(s)
end


"Run pipeline from command line."
function runfromshell(parsed_args)
    runpipeline(;
        i1=parsed_args["i1"],
        d1=parsed_args["d1"],
        i2=parsed_args["i2"],
        d2=parsed_args["d2"],
        out_dir=parsed_args["out_dir"],
        n_strucs=parsed_args["n_strucs"],
        tolerance_weight=parsed_args["tolerance_weight"],
        other_ratio=parsed_args["other_ratio"],
        extra_pdbs=parsed_args["extra_pdbs"],
        mod_path=parsed_args["mod_path"],
        n_mods=parsed_args["n_mods"]
    )
end


"Run auto-parameterisation from command line."
function paramfromshell(parsed_args)
    parampipeline(;
        i1=parsed_args["i1"],
        d1=parsed_args["d1"],
        i2=parsed_args["i2"],
        d2=parsed_args["d2"],
        out_dir=parsed_args["out_dir"],
        n_strucs=parsed_args["n_strucs"],
        other_ratio=parsed_args["other_ratio"],
        frac_between=parsed_args["fraction"],
        tmscore_path=parsed_args["tmscore"]
    )
end


"Wrapper function to run the whole pipeline."
function runpipeline{T <: AbstractString}(;
                    i1::Union{AbstractString, Void}=nothing,
                    d1::Union{AbstractString, Void}=nothing,
                    i2::Union{AbstractString, Void}=nothing,
                    d2::Union{AbstractString, Void}=nothing,
                    out_dir::AbstractString=defaults["out_dir"],
                    n_strucs::Integer=defaults["n_strucs"],
                    tolerance_weight::Real=defaults["tolerance_weight"],
                    other_ratio::Real=defaults["other_ratio"],
                    extra_pdbs::Array{T,1}=String[],
                    mod_path::Union{AbstractString, Void}=nothing,
                    n_mods::Integer=0)
    @assert i1 != nothing && d1 != nothing "Arguments i1 and d1 required"

    # Print inputs
    println()
    println("-- ExProSE --")
    println()
    println("Arguments:")
    println("  i1               - ", i1)
    println("  d1               - ", d1)
    println("  i2               - ", i2)
    println("  d2               - ", d2)
    println("  out_dir          - ", out_dir)
    println("  n_strucs         - ", n_strucs)
    println("  tolerance_weight - ", tolerance_weight)
    println("  other_ratio      - ", other_ratio)
    if length(extra_pdbs) > 0
        println("  extra_pdbs       - ", join(extra_pdbs, ", "))
    else
        println("  extra_pdbs       - nothing")
    end
    println("  mod_path         - ", mod_path)
    println("  n_mods           - ", n_mods)
    println()

    # Make output directories
    inner_dirs = ["pdbs", "pymol", "flucs", "pcs", "scores", "rmsds"]
    append!(inner_dirs, ["pdbs_mod_$i" for i in 1:n_mods])
    append!(inner_dirs, ["mod_$i" for i in 1:n_mods])
    makedirectories(out_dir, inner_dirs)

    # Decide whether to calculate the constraints from one or two structures
    if i2 != nothing && d2 != nothing
        println("Proceeding with two structures")
        constraints_com, constraints_one, constraints_two = interactions(i1, d1, i2, d2, other_ratio=other_ratio, tolerance_weight=tolerance_weight)
        ensemble_com = generateensemble(constraints_com, n_strucs)
        # Perturbation carried out on holo structure
        ensemble_mods = perturbensemble(constraints_two.atoms, constraints_com, n_strucs, mod_path, n_mods)
        runanalysis(out_dir, ensemble_com, constraints_one, constraints_two, extra_pdbs=extra_pdbs, ensemble_mods=ensemble_mods)
    else
        println("Proceeding with one structure")
        constraints = interactions(i1, d1, other_ratio=other_ratio, tolerance_weight=tolerance_weight)
        ensemble = generateensemble(constraints, n_strucs)
        ensemble_mods = perturbensemble(constraints.atoms, constraints, n_strucs, mod_path, n_mods)
        runanalysis(out_dir, ensemble, constraints, extra_pdbs=extra_pdbs, ensemble_mods=ensemble_mods)
    end
    println("Done")
end


"Make the outer results directory and any inner directories from a list."
function makedirectories(out_dir::AbstractString, inner_dirs::Array{String,1}=String[])
    if isdir(out_dir)
        println("Output directory \"$out_dir\" already exists - beware of overwriting files!")
    else
        println("Making output directory \"$out_dir\"")
        mkdir(out_dir)
    end
    for inner_dir in inner_dirs
        new_dir = "$(out_dir)/$inner_dir"
        if !isdir(new_dir)
            mkdir(new_dir)
        end
    end
end


"Run the analysis pipeline and write output files."
function runanalysis{T <: AbstractString}(
                    out_dir::AbstractString,
                    ensemble::ModelledEnsemble,
                    constraints::Constraints;
                    extra_pdbs::Array{T,1}=String[],
                    ensemble_mods::Array{ModelledEnsemble,1}=ModelledEnsemble[],
                    out_prefix::AbstractString=defaults["out_prefix"])
    # Align ensemble
    ens_al = selfalignensemble!(ensemble)
    alignatoms!(constraints.atoms, ens_al)
    ens_av = centroid(ensemble)

    # Do PCA and write projections
    pca = PCA(ensemble)
    pcs_ref = projectstructure(constraints.atoms, pca)
    writeprojections("$out_dir/pcs/pcs.dat", pca.pcs)
    writeprojections("$out_dir/pcs/pcs_ref.dat", pcs_ref)
    writefloatarray("$out_dir/pcs/evals_spread.txt", pca.evals / sum(pca.evals))

    # Write PDBs, SPE scores and RMSDs
    writeensemble("$out_dir/pdbs/$out_prefix.pdb", ensemble)
    writeensemblescores("$out_dir/scores/scores.txt", ensemble)
    writefloatarray("$out_dir/rmsds/rmsds.txt", ensemblermsds(ensemble, constraints.atoms))

    # Plot the PCs
    if length(extra_pdbs) > 0
        pcs_extra = projectpdbfiles("$out_dir/pcs/pcs_ex_", extra_pdbs, constraints.atoms, ens_al, pca)
    else
        pcs_extra = Float64[]
    end
    plotpcs("$out_dir/pcs/pc", pca.pcs, pcs_ref_one=pcs_ref, pcs_extra=pcs_extra)

    # Write and plot ensemble fluctuations
    flucs = fluctuations(ensemble)
    writefloatarray("$out_dir/flucs/flucs.txt", flucs)
    plotfluctuations("$out_dir/flucs/flucs.png", flucs)

    # Deal with perturbed ensembles
    av_rmsds = Float64[]
    for (i, ensemble_mod) in enumerate(ensemble_mods)
        alignensemble!(ensemble_mod, ens_al)
        push!(av_rmsds, rmsd(ens_av, centroid(ensemble_mod), constraints.atoms))
        pcs_mod = projectensemble(ensemble_mod, pca)
        writeprojections("$out_dir/pcs/pcs_mod_$i.dat", pcs_mod)
        writeensemble("$out_dir/pdbs_mod_$i/$out_prefix.pdb", ensemble_mod)
        writeensemblescores("$out_dir/scores/scores_mod_$i.txt", ensemble_mod)
        plotpcs("$out_dir/mod_$i/pc", pca.pcs, pcs_ref_one=pcs_ref, pcs_mod=pcs_mod, pcs_extra=pcs_extra)
        writefloatarray("$out_dir/mod_$i/rmsds.txt", ensemblermsds(ensemble_mod, constraints.atoms))
        selfalignensemble!(ensemble_mod)
        flucs_mod = fluctuations(ensemble_mod)
        writefloatarray("$out_dir/flucs/flucs_mod_$i.txt", flucs_mod)
        writefloatarray("$out_dir/flucs/fluc_rats_mod_$i.txt", flucs_mod ./ flucs)
        plotfluctuations("$out_dir/mod_$i/flucs.png", flucs, flucs_mod=flucs_mod)
    end

    # Write PyMol scripts to view PCs
    writepcviews("$out_dir/pymol/pc_", pca, ensemble.atoms, 5)

    # Write input atoms back out
    writepdb("$out_dir/input.pdb", constraints.atoms)

    # Write perturbation values
    if length(ensemble_mods) > 0
        writefloatarray("$out_dir/perturbations.txt", av_rmsds)
    end
end


function runanalysis{T <: AbstractString}(
                    out_dir::AbstractString,
                    ensemble_com::ModelledEnsemble,
                    constraints_one::Constraints,
                    constraints_two::Constraints;
                    extra_pdbs::Array{T,1}=String[],
                    ensemble_mods::Array{ModelledEnsemble,1}=ModelledEnsemble[],
                    out_prefix::AbstractString=defaults["out_prefix"])
    # Align ensemble
    ens_al = selfalignensemble!(ensemble_com)
    alignatoms!(constraints_one.atoms, ens_al)
    alignatoms!(constraints_two.atoms, ens_al)
    ens_av = centroid(ensemble_com)

    # Do PCA and write projections
    pca = PCA(ensemble_com)
    pcs_ref_one = projectstructure(constraints_one.atoms, pca)
    pcs_ref_two = projectstructure(constraints_two.atoms, pca)
    writeprojections("$out_dir/pcs/pcs_com.dat", pca.pcs)
    writeprojections("$out_dir/pcs/pcs_ref_one.dat", pcs_ref_one)
    writeprojections("$out_dir/pcs/pcs_ref_two.dat", pcs_ref_two)
    writeintarray("$out_dir/pcs/pcs_ord.txt", findimportantpcs(pcs_ref_one, pcs_ref_two))
    writefloatarray("$out_dir/pcs/evals_spread.txt", pca.evals / sum(pca.evals))

    # Write PDBs, SPE scores and RMSDs
    writeensemble("$out_dir/pdbs/$out_prefix.pdb", ensemble_com)
    writeensemblescores("$out_dir/scores/scores.txt", ensemble_com)
    writefloatarray("$out_dir/rmsds/rmsds_ref_one.txt", ensemblermsds(ensemble_com, constraints_one.atoms))
    writefloatarray("$out_dir/rmsds/rmsds_ref_two.txt", ensemblermsds(ensemble_com, constraints_two.atoms))

    # Plot the PCs
    if length(extra_pdbs) > 0
        pcs_extra = projectpdbfiles("$out_dir/pcs/pcs_ex_", extra_pdbs, constraints_one.atoms, ens_al, pca)
    else
        pcs_extra = Float64[]
    end
    plotpcs("$out_dir/pcs/pc", pca.pcs, pcs_ref_one=pcs_ref_one, pcs_ref_two=pcs_ref_two, pcs_extra=pcs_extra)

    # Write and plot ensemble fluctuations
    flucs = fluctuations(ensemble_com)
    writefloatarray("$out_dir/flucs/flucs.txt", flucs)
    plotfluctuations("$out_dir/flucs/flucs.png", flucs)

    # Deal with perturbed ensembles
    av_rmsds = Float64[]
    for (i, ensemble_mod) in enumerate(ensemble_mods)
        alignensemble!(ensemble_mod, ens_al)
        push!(av_rmsds, rmsd(ens_av, centroid(ensemble_mod), constraints_one.atoms))
        pcs_mod = projectensemble(ensemble_mod, pca)
        writeprojections("$out_dir/pcs/pcs_mod_$i.dat", pcs_mod)
        writeensemble("$out_dir/pdbs_mod_$i/$out_prefix.pdb", ensemble_mod)
        writeensemblescores("$out_dir/scores/scores_mod_$i.txt", ensemble_mod)
        plotpcs("$out_dir/mod_$i/pc", pca.pcs, pcs_ref_one=pcs_ref_one, pcs_ref_two=pcs_ref_two, pcs_mod=pcs_mod, pcs_extra=pcs_extra)
        writefloatarray("$out_dir/mod_$i/rmsds_ref_one.txt", ensemblermsds(ensemble_mod, constraints_one.atoms))
        writefloatarray("$out_dir/mod_$i/rmsds_ref_two.txt", ensemblermsds(ensemble_mod, constraints_two.atoms))
        selfalignensemble!(ensemble_mod)
        flucs_mod = fluctuations(ensemble_mod)
        writefloatarray("$out_dir/flucs/flucs_mod_$i.txt", flucs_mod)
        writefloatarray("$out_dir/flucs/fluc_rats_mod_$i.txt", flucs_mod ./ flucs)
        plotfluctuations("$out_dir/mod_$i/flucs.png", flucs, flucs_mod=flucs_mod)
    end

    # Write PyMol scripts to view PCs
    writepcviews("$out_dir/pymol/pc_", pca, ensemble_com.atoms, 5)

    # Write input atoms back out
    writepdb("$out_dir/input_1.pdb", constraints_one.atoms)
    writepdb("$out_dir/input_2.pdb", constraints_two.atoms)

    # Write perturbation values
    if length(ensemble_mods) > 0
        writefloatarray("$out_dir/perturbations.txt", av_rmsds)
    end
end


"Run the auto-parameterisation pipeline."
function parampipeline(;
                    i1::Union{AbstractString, Void}=nothing,
                    d1::Union{AbstractString, Void}=nothing,
                    i2::Union{AbstractString, Void}=nothing,
                    d2::Union{AbstractString, Void}=nothing,
                    out_dir::AbstractString=defaults["out_dir_param"],
                    n_strucs::Integer=defaults["n_strucs_param"],
                    other_ratio::Real=defaults["other_ratio"],
                    frac_between::Real=defaults["frac_between"],
                    tmscore_path::AbstractString=defaults["tmscore_path"],
                    out_prefix::AbstractString=defaults["out_prefix"],
                    param_tw_start::Real=defaults["param_tw_start"],
                    param_tw_increment::Real=defaults["param_tw_increment"])
    @assert i1 != nothing && d1 != nothing && i2 != nothing && d2 != nothing "Arguments i1, d1, i2 and d2 required"
    @assert 0.0 <= frac_between <= 1.0 "frac_between cannot be less than 0.0 or more than 1.0"
    @assert param_tw_start >= 0.0 "param_tw_start cannot be negative"
    @assert param_tw_increment > 0.0 "param_tw_increment must be positive"
    @assert tmscorepathvalid(tmscore_path) "Not a valid TMscore path: \"$tmscore_path\""

    # Print inputs
    println()
    println("-- ExProSE auto-parameterisation --")
    println()
    println("Arguments:")
    println("  i1               - ", i1)
    println("  d1               - ", d1)
    println("  i2               - ", i2)
    println("  d2               - ", d2)
    println("  out_dir          - ", out_dir)
    println("  n_strucs         - ", n_strucs)
    println("  other_ratio      - ", other_ratio)
    println("  frac_between     - ", frac_between)
    println("  tmscore          - ", tmscore_path)
    println()
    println("Beginning auto-parameterisation")

    makedirectories(out_dir)
    tw = param_tw_start
    found = false
    fracs = Float64[]
    while !found
        if tw < 0.0
            println("Tolerance weight cannot be reduced any further")
            println("Suggested tolerance weight is 0.0")
            break
        end
        println("Running with tolerance weight of ", tw)
        constraints_com, constraints_one, constraints_two = interactions(i1, d1, i2, d2, other_ratio=other_ratio, tolerance_weight=tw)
        ensemble_com = generateensemble(constraints_com, n_strucs)
        selfalignensemble!(ensemble_com)
        lab = replace(string(tw), ".", "_") # Replace dots with underscores for directory name
        if !isdir("$out_dir/tw_$lab")
            mkdir("$out_dir/tw_$lab")
            mkdir("$out_dir/tw_$lab/pdbs")
        end
        ens_prefix = "$out_dir/tw_$lab/pdbs/$out_prefix.pdb"
        writeensemble(ens_prefix, ensemble_com)
        frac = fractionbetweeninputs(ens_prefix, length(ensemble_com.strucs), i1, i2; command=tmscore_path)
        println("Fraction of structures between inputs is ", frac)
        push!(fracs, frac)
        if frac >= frac_between
            println("Fraction is at least threshold of ", frac_between)
            println("Suggested tolerance weight is ", tw)
            found = true
        else
            println("Fraction is below threshold of ", frac_between)
            println("Reducing tolerance weight and running again")
        end
        # This causes weird float things to happen
        tw -= param_tw_increment
    end
    # Should write tw and fraction
    writefloatarray("$out_dir/fractions.txt", fracs; dec_places=2)
    println("Done")
end
