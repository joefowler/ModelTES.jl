using ModelTES, PyPlot

"""`noise_power_demo()`

Plot a noise power spectral density for five TES model systems.
For one, plot the 4 components from the Irwin-Hilton model."""

function noise_power_demo(SI_amp=5e-22)
    models = Dict(
        :LowE => ModelTES.lowEpix(),
        :HighE => ModelTES.highEpix(),
        :Holmes12 => ModelTES.pholmes(12e-9),
        :Holmes24 => ModelTES.pholmes(24e-9),
        :Holmes48 => ModelTES.pholmes(48e-9))
    colors=Dict(:HighE=>"r", :LowE=>"gold",
        :Holmes12=>"purple",:Holmes24=>"blue",:Holmes48=>"green");

    clf()
    freq = logspace(1, 6, 50);
    for m in (:LowE, :HighE, :Holmes12, :Holmes24, :Holmes48)
        model = IrwinHiltonTES(models[m])
        psd,N1,N2,N3,N4 = ModelTES.noisePSD(model, freq, SI_amp)
        loglog(freq,psd, color=colors[m], label = string(m))
        if m==:HighE
            loglog(freq, N1, "--", color=colors[m], label="HighE TES")
            loglog(freq, N2, ":", color=colors[m], label="HighE Load")
            loglog(freq, N3, "-.", color=colors[m], label="HighE TFN")
            loglog(freq, N4, "--", color="orange", label="HighE ampl")
        end
    end
    xlabel("Frequency (Hz)")
    ylabel("Noise power (A\$^2\$/Hz)")
    title("Noise Power Spectrum for several TES models")
    ylim([.3*SI_amp,3e-20+SI_amp])
    legend(loc="best")
    grid(true)
end
