module DataScripts

using PyPlot
using DataFrames
using SFGTools

export plot

const DATA_PATH = "/Users/lackner/Documents/DataSFG"

# f = figure()
# close(f)

function plot(idx=1::Int64; pixel=false, sfwl=false)
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

end #module
