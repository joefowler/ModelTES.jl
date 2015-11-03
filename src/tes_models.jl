using PyPlot, ModelTES

"""Return a BiasedTES object containing parameters for a HOLMES notional TES"""
function pholmes(L)
    n = 3.25; Tc = 0.1; Tbath = 0.07
    k = 2.33e-8; C = 0.5e-12
    R0 = 2e-3; Rl = 0.3e-3; Rn = 10e-3; Rpara=0.0
    alpha = 200.0; beta = 2.0
    model = ModelTES.ShankRIT(alpha, beta, n, Tc, Tbath, k, R0, Rn);
    tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
    bt = BiasedTES(tes_param, R0)
end


"""Return a BiasedTES object containing parameters for Kelsey's low-E TES model"""
function lowEpix()
    Tc=0.094; Tbath = 0.075
    R0 = 1.64e-3; Rl = 0.35e-3; Rn = 8.2e-3; Rpara=0.0
    n=3.15; k = 1.48e-9; C = 0.1e-12
    alpha = 175.0; beta = 2.28; L=50e-9
    model = ModelTES.ShankRIT(alpha, beta, n, Tc, Tbath, k, R0, Rn);
    tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
    bt = BiasedTES(tes_param, R0)
end


"""Return a BiasedTES object containing parameters for Kelsey's high-E TES model"""
function highEpix()
    Tc=0.107; Tbath = 0.065
    R0 = 1.55e-3; Rl = 0.35e-3; Rn = 10.3e-3; Rpara=0.0
    n=3.3; k = 3.50e-9; C = 0.72e-12
    alpha = 92.0; beta = 1.4; L=300e-9
    model = ModelTES.ShankRIT(alpha, beta, n, Tc, Tbath, k, R0, Rn);
    tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
    bt = BiasedTES(tes_param, R0)
end

"Plot a noise power spectral density for 5 TES model systems.
For one, plot the 4 components from the Irwin-Hiltom model."
function noise_power_demo(SI_amp=5e-22)
    models = Dict(
        :LowE=>lowEpix(),
        :HighE=>highEpix(),
        :Holmes12=>pholmes(12e-9),
        :Holmes24=>pholmes(24e-9),
        :Holmes48=>pholmes(48e-9))
    colors=Dict(:HighE=>"r", :LowE=>"gold",
        :Holmes12=>"purple",:Holmes24=>"blue",:Holmes48=>"green");

    clf()
    freq = logspace(1, 6, 50);
    for m in (:LowE, :HighE, :Holmes12, :Holmes24, :Holmes48)
        model = IrwinHiltonTES(models[m])
        psd,N1,N2,N3,N4 = noise(model, freq, SI_amp)
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
