# Run TMscore from command line - see https://zhanglab.ccmb.med.umich.edu/TM-score


export
    tmscore,
    tmscorepathvalid,
    fractionbetweeninputs


"""
TM-score of model and reference, given as filepaths.
Uses external command line program - see https://zhanglab.ccmb.med.umich.edu/TM-score.
"""
function tmscore(model::AbstractString,
                    reference::AbstractString;
                    command::AbstractString=defaults["tmscore_path"])
    @assert isfile(model) "Model filepath not valid: \"$model\""
    @assert isfile(reference) "Reference filepath not valid: \"$reference\""
    tmscore_output = readstring(pipeline(`TMscore $model $reference`))
    # This is easier using grep and sed but they are system-specific
    reg_match = match(r"\nTM-score    = \d\.\d\d\d\d", tmscore_output)
    @assert reg_match != nothing "TMscore failed or the result could not be read"
    return float(reg_match.match[16:end])
end


"Check a command runs TMscore."
function tmscorepathvalid(command::AbstractString)
    try
        readstring(pipeline(`$command`))
    catch
        return false
    end
    return true
end


"""
Fraction of structures in an ensemble that are closer to both input structures
than the input structures are to each other.
"""
function fractionbetweeninputs(ens_prefix::AbstractString,
                    n_strucs::Integer,
                    i1::AbstractString,
                    i2::AbstractString;
                    command::AbstractString=defaults["tmscore_path"])
    tmscores_ens_i1 = [tmscore("$ens_prefix.$i", i1; command=command) for i in 1:n_strucs]
    tmscores_ens_i2 = [tmscore("$ens_prefix.$i", i2; command=command) for i in 1:n_strucs]
    tmscore_i1_i2 = tmscore(i1, i2; command=command)
    tmscore_i2_i1 = tmscore(i2, i1; command=command)
    counter = 0
    for i in 1:n_strucs
        if tmscores_ens_i1[i] > tmscore_i2_i1 && tmscores_ens_i2[i] > tmscore_i1_i2
            counter += 1
        end
    end
    return counter / n_strucs
end
