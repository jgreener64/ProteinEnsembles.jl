# Plot results


export
    plotpcs,
    plotfluctuations


using Gadfly
using Colors
import Cairo


"Dictionary of default colours for plotting."
const default_colours = Dict{String, String}(
    "exprose"=> "#3498DB", # Blue
    "input_1"=> "#E74C3C", # Red
    "input_2"=> "#2ECC71", # Green
    "mod"    => "#FF9900", # Orange
    "extra"  => "#000000", # Black
)


"`Theme` for plots."
plottheme() = Theme(
    background_color=colorant"white",
    panel_stroke=colorant"black",
    major_label_color=colorant"black",
    minor_label_color=colorant"black",
    major_label_font_size=20pt,
    minor_label_font_size=16pt,
    grid_line_width=0mm,
)


"`Theme` for plot points."
pointtheme(colour::AbstractString) = Theme(
    default_color=parse(Colorant, colour),
    point_size=1.5mm,
    highlight_width=0mm,
)


"Plot the projections of an ensemble onto the principal components of the motion."
function plotpcs(out_prefix::AbstractString,
                    pcs::Array{Float64};
                    pcs_ref_one::Array{Float64}=Float64[],
                    pcs_ref_two::Array{Float64}=Float64[],
                    pcs_mod::Array{Float64}=Float64[],
                    pcs_extra::Array{Float64}=Float64[],
                    comps_to_plot::Array{Tuple{Int,Int},1}=[(1,2), (1,3), (2,3)])
    checkfilepath(out_prefix)
    for (i, j) in comps_to_plot
        out_filepath = "$(out_prefix)_$(i)_$j.png"
        @assert length(pcs) > 0 "No projections in array"
        @assert max(i, j) <= size(pcs, 1) "Chose a PC number higher than those available"
        layers = layer(
            x=pcs[i,:],
            y=pcs[j,:],
            Geom.point,
            pointtheme(default_colours["exprose"]),
        )
        pcsaddlayer!(layers, pcs_ref_one, (i, j), default_colours["input_1"])
        pcsaddlayer!(layers, pcs_ref_two, (i, j), default_colours["input_2"])
        pcsaddlayer!(layers, pcs_mod, (i, j), default_colours["mod"])
        pcsaddlayer!(layers, pcs_extra, (i, j), default_colours["extra"])
        p = plot(
            layers...,
            Guide.xlabel("Projections onto PC $i / Å"),
            Guide.ylabel("Projections onto PC $j / Å"),
            plottheme,
        )
        draw(PNG(out_filepath, 20cm, 20cm), p)
    end
    println("Plotted principal components to file(s) \"$(out_prefix)_x_y.png\"")
end


"Adds a layer of points to an existing layer for PC plotting."
function pcsaddlayer!(layers::Array{Layer,1}, pcs::Array{Float64}, comps::Tuple{Int,Int}, colour::AbstractString)
    if length(pcs) > 0
        i, j = comps
        @assert max(i, j) <= size(pcs, 1) "Chose a PC number higher than those available"
        push!(layers, layer(
            x=pcs[i,:],
            y=pcs[j,:],
            Geom.point,
            pointtheme(colour),
        )[1])
    end
end


"""
Plot root mean square fluctuations of an aligned ensemble.
Arguments are the output filepath and the array of fluctuations.
"""
function plotfluctuations(out_filepath::AbstractString,
                    flucs::Array{Float64,1};
                    flucs_mod::Array{Float64,1}=Float64[])
    checkfilepath(out_filepath)
    n_res = length(flucs)
    @assert n_res > 0 "Fluctuation array is empty"
    layers = layer(
        x=collect(1:n_res),
        y=flucs,
        Geom.line,
        pointtheme(default_colours["exprose"]),
    )
    if length(flucs_mod) > 0
        @assert length(flucs_mod) == n_res "Number of fluctuations differ"
        push!(layers, layer(
            x=collect(1:n_res),
            y=flucs_mod,
            Geom.line,
            pointtheme(default_colours["mod"]),
        )[1])
    end
    p = plot(
        layers...,
        Guide.xlabel("Residue index"),
        Guide.ylabel("RMSF / Å"),
        plottheme,
    )
    draw(PNG(out_filepath, 20cm, 15cm), p)
    println("Plotted RMSFs to file \"$out_filepath\"")
end
