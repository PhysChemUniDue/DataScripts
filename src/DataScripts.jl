module DataScripts

using LsqFit
using Plots
using DataFrames
using SFGTools

export plotrecent, pump_probe_process, cal_spectrometer

const DATA_PATH = "/Users/lackner/Documents/DataSFG"

# f = figure()
# close(f)

function plotrecent(idx=1; pixel=false, sfwl=false, dograb=false)
    dograb && grab(DATA_PATH)
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

    data = average.(data)

    # Empty Plot
    plot()
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
        plot!(x, data[i].s, label=name)
        xlabel!(xlabelstr)
        ylabel!("Counts per Second")
    end
    gui()
    return nothing
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
    data = SFGTools.average.(data)

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


function cal_spectrometer(idarray)

    data = SFGTools.load_spectra(idarray)
    data = SFGTools.average.(data)
    names = get_attribute(data, "name")
    splitted_names = split.(names, "_")

    λcal = zeros(length(data))
    locs = fill("", length(data))
    for i in eachindex(splitted_names)
        λcal[i] = splitted_names[i][end-1] |> parse
        locs[i] = splitted_names[i][end]
    end

    # Check if the data is suitable
    for i = 1:length(data)
        length(data[i].s) != 512 && error("All spectra need to have 512 data points")
        !any(["center", "right", "left"] .== locs[i]) && error("Filenames have to have the form \"*_334.15_left\" etc. (see help)")
    end

    # Remove everything but the peaks
    for i in eachindex(data)
        if locs[i] == "center"
            data[i].s[[1:236..., 276:end...]] = NaN
        elseif locs[i] == "left"
            data[i].s[20:end] = NaN
        elseif locs[i] == "right"
            data[i].s[1:end-20] = NaN
        end
    end

    peakpositions = fill(0, length(data))
    for i in 1:length(data)
        peakpositions[i] = findmax(data[i].s)[2]
    end
    peakpositions = reshape(peakpositions, (3, length(data) ÷ 3))'
    peakpositions = hcat(peakpositions[:,2], peakpositions[:,1], peakpositions[:,3])

    λ_spectrometer = Float64[]
    for (i, d) in enumerate(data)
        λ = get_attribute(d, "spectrometer_wavelength")
        push!(λ_spectrometer, λ)
    end
    λ_spectrometer = reshape(λ_spectrometer, (3, 7))'
    λ_spectrometer = hcat(λ_spectrometer[:,2], λ_spectrometer[:,1], λ_spectrometer[:,3])

    Δp = abs.(diff(peakpositions, 2))
    Δλ = abs.(diff(λ_spectrometer, 2))
    nmperpixel = Δλ ./ Δp

    mean_nmperpixel = mean(nmperpixel, 2)
    mean_nmperpixel = squeeze(mean_nmperpixel, 2)

    model(x, p) = p[1] * x.^2 + p[2] * x + p[3]
    fit = curve_fit(model, λ_spectrometer[:,2], mean_nmperpixel, [-3e-8, 5e-6, 0.02])
    xfit = linspace(330, 600, 100)
    yfit = model(xfit, fit.param)

    p1 = scatter(λ_spectrometer[:,2], mean_nmperpixel)
    plot!(xfit, yfit)
    xlabel!("Central Spectrometer Wavelength in nm")
    ylabel!("Wavelength Density [nm/pixel]")

    # the real center is in between pixels 256 and 257!
    Δp0 = peakpositions[:,2] - 256.5

    offset_fit = curve_fit(model, λ_spectrometer[:,2], Δp0, [-3e-4, 0.25, -32.0])
    offset_yfit = model(xfit, offset_fit.param)

    p2 = scatter(λ_spectrometer[:,2], Δp0)
    plot!(xfit, offset_yfit)
    xlabel!("Central Spectrometer Wavelength in nm")
    ylabel!("Pixel Offset")

    plot(p1, p2)
    gui()

    info("Add the following code to the `get_wavelength` function at the appropriate location in the SFGTools package:")

    println("
    # Calibration parameters were determined on $(now())
    # Pixel offset
    Δp = $(offset_fit.param[1]) * λ0^2 + $(offset_fit.param[2]) * λ0 - $(offset_fit.param[3])
    # Wavelenths per pixel
    pλ = $(fit.param[1]) * λ0^2 + $(fit.param[2]) * λ0 + $(fit.param[3])
    ")

    info("The right location is inside the if-elseif-else statement. Change the last else to:")

    println("elseif date < DateTime(\"$(now())\")")

    info("and add a new final else statement with the above code that contains the calibration parameters.")
    return nothing
end

end #module
