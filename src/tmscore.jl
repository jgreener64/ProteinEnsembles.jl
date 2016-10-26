# Run TMscore from command line - see http://zhanglab.ccmb.med.umich.edu/TM-score


export
    tmscore,
    tmscorepathvalid


""
function tmscore()

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
