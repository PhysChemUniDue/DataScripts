module DataScripts

using Plots
using DataFrames
import SFGTools

export plotrecent, pump_probe_process

const DATA_PATH = "/Users/lackner/Documents/DataSFG"

# f = figure()
# close(f)

function plotrecent(idx=1::Int64; pixel=false, sfwl=false)
    grab(DATA_PATH)
    df = list_spectra()
    sort!(df, cols=:id, rev=true)

    if maximum(idx) < 9999
        # Treat as DataFrame Index
        data = load_spectra(df[:id][idx])
    else
        # Treat as Spectrum ID
        data = load_spectra(idx)
    end

    if length(idx) == 1
        # Make data iterable
        data = [data]
    end

    figure()
    for i = 1:length(data)
        if pixel
            x_binning = get_attribute(data[i], "x_binning")
            npixel = 512 / x_binning
            x = 1:npixel
            xlabelstr = "Pixel Number"
        elseif sfwl
            x = get_wavelength(data[i])
            xlabelstr = "Wavelength in nm"
        else
            x = get_ir_wavelength(data[i])
            xlabelstr = "IR Wavelength in nm"
        end

        name = get_attribute(data[i], "name")
        plot(x, mean(data[i]), label=name)
        xlabel(xlabelstr)
        ylabel("Counts per Second")
    end
    legend(loc="lower center", bbox_to_anchor=(0.5, 1.0),
          fancybox=false, shadow=false, ncol=3)
end

"""
Process pump probe spectra. The data is splitted in the middle of the camera.

Keyword Arguments:

* `eventwidth`: Max. width of events that are being removed (default 8)
* `savepath`:   Path were the processed data is stored (default is the folder where the raw data is loaded from)
* `wlrange`:    Figures are cropped to these xlims (defaults to full range of available values) eg.: `wlrange=(3300, 3600)`
"""
function pump_probe_process(idarray; 
    eventwidth=8, 
    savepath="default",
    wlrange="default")

    print("Loading Data...")
    data = SFGTools.load_spectra(idarray)
    SFGTools.rm_events!.(data)
    SFGTools.mean!.(data)

    if savepath == "default"
        savepath = splitdir(splitdir(SFGTools.getdir(data[1].id))[1])[1]
        savepath = joinpath(savepath, "processed")
    end
    !isdir(savepath) && mkdir(savepath)

    dltime = SFGTools.get_pump_delay.(data)
    λ = SFGTools.get_ir_wavelength(data[1])
    ν = SFGTools.get_ir_wavenumber(data[1])
    name = SFGTools.get_attribute(data[1], "name")

    # Sort with respect to delay times
    p = sortperm(dltime)
    dltime = dltime[p]
    data = data[p]
    println("\tDone.")

    # Check if all delays are unique
    length(unique(dltime)) != length(dltime) && warn("Not all delay times are unique!")

    # Get Wavelength range
    wlrange == "default" && (wlrange = (minimum(λ), maximum(λ)))

    # Make seperate spectra for upper and lower part of CCD camera for all time
    # delays.
    print("Generating Spectra...")
    npixels = size(data[1].s)
    spectra = zeros(npixels[1],2)
    split = npixels[2] ÷ 2
    for d in data
        # Upper part
        spectra[:,1] += squeeze(sum(d.s[:,1:split], 2), 2)
        # Lower part
        spectra[:,2] += squeeze(sum(d.s[:,split+1:end], 2), 2)
    end
    # Normalize
    spectra /= length(data)

    plotarray = Plots.Plot[]
    p1 = Plots.plot(λ, spectra, labels=["A", "B"])
    Plots.title!(name * "/nSpectra Total")
    Plots.ylabel!("Counts [1/s]")
    Plots.xlabel!("IR Wavelength [nm]")
    Plots.xlims!(wlrange)
    push!(plotarray, p1)
    Plots.savefig(joinpath(savepath, "spectra.png"))
    println("\tDone.")

    # Make Heatmap
    print("Generating Heatmaps...")
    hmap = zeros(length(data), npixels[1], 2)
    for i in eachindex(data)
        hmap[i,:,1] = data[i].s[:,1:split]
        hmap[i,:,2] = data[i].s[:,split+1:end]
    end
    p2 = Plots.heatmap(λ, dltime, hmap[:,:,1], clim=(minimum(hmap), maximum(hmap)))
    Plots.title!(name * "\nA")
    Plots.xlabel!("IR Wavelength [nm]")
    Plots.ylabel!("Pump Delay [ps]")
    Plots.xlims!(wlrange)
    push!(plotarray, p2)
    Plots.savefig(joinpath(savepath, "heatmap_A.png"))
    
    p3 = heatmap(λ, dltime, hmap[:,:,2], clim=(minimum(hmap), maximum(hmap)))
    Plots.title!(name * "\nB")
    Plots.xlabel!("IR Wavelength [nm]")
    Plots.ylabel!("Pump Delay [ps]")
    Plots.xlims!(wlrange)
    push!(plotarray, p3)
    Plots.savefig(joinpath(savepath, "heatmap_B.png"))
    println("\tDone.")

    # Animated Spectra
    print("Animating Spectra")
    labels=("A", "B")
    anim = Plots.@animate for i = 1:size(hmap, 1)
        print(".")
        
        p1 = Plots.plot(λ, hmap[i,:,:], labels=labels, ylim=(minimum(hmap), maximum(hmap)))
        Plots.title!("Delay = $(round(dltime[i], 1)) ps")
        Plots.ylabel!("Counts [1/s]")
        Plots.xlabel!("IR Wavelength [nm]")
        Plots.xlims!(wlrange)

        p2 = Plots.heatmap(λ, dltime, hmap[:,:,1], clim=(minimum(hmap), maximum(hmap)))
        Plots.title!(name * "\nA")
        Plots.xlabel!("IR Wavelength [nm]")
        Plots.ylabel!("Pump Delay [ps]")
        Plots.scatter!([wlrange[1], wlrange[end]], [dltime[i], dltime[i]], markersize=10, legend=false, color=:white)
        Plots.xlims!(wlrange)

        p3 = heatmap(λ, dltime, hmap[:,:,2], clim=(minimum(hmap), maximum(hmap)))
        Plots.title!(name * "\nB")
        Plots.xlabel!("IR Wavelength [nm]")
        Plots.ylabel!("Pump Delay [ps]")
        Plots.scatter!([wlrange[1], wlrange[end]], [dltime[i], dltime[i]], markersize=10, legend=false, color=:white)
        Plots.xlims!(wlrange)

        Plots.plot(p1, p2, p3, layout=(3,1), size=(600,1200))
    end
    Plots.gif(anim, joinpath(savepath, "spectra.gif"), fps=2)
    println("\tDone.")

    return spectra, plotarray
end

end #module
