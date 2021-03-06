using ModelTES, ARMA
using Base.Test

# Create a high-E TES design
stdTES = ModelTES.highEpix()
# iv_pt = ModelTES.iv_point(stdTES.p, stdTES.V)
Vs_in = stdTES.V*collect(0:0.1:10)
Is, Ts, Rs, Vs_out = ModelTES.iv_curve(stdTES.p, Vs_in)
derivs = [ModelTES.dT_and_dI_iv_point(stdTES.p, Is[i], Ts[i], Rs[i], Vs_in[i]) for i in eachindex(Vs_in)]
@test maximum(map(maximum,abs.(derivs)))<1e-10

# @show iv_curve
# Create a Biased TES from the 48 nanohentry Holmes paramters with 0.2*Rn resistance
tes = ModelTES.pholmes(48e-9, 0.20)
# Create a Biased TES from the 48 nanohentry Holmes paramters with 0.4*Rn resistance
tes2 = ModelTES.pholmes(48e-9, 0.40)
# Integrate a pulse with 12000 samples, 1e-7 second spacing, 1000 eV energy, 2000 presamples
out = rk8(12000,1e-7, tes, 1000, 2000);
out_pulse = pulse(12000,1e-7, tes, 1000, 2000);
function worst_relative_error(a,b)
    @assert(all(times(a).==times(b)))
    eI=maximum(abs.(2*(a.I-b.I)./(a.I+b.I)))
    eT=maximum(abs.(2*(a.T-b.T)./(a.T+b.T)))
    eR=maximum(abs.(2*(a.R-b.R)./(a.R+b.R)))
    max(eI,eT,eR)
end
@test worst_relative_error(out,out_pulse)<1e-5
out_temp = deepcopy(out)
out_temp.I[1]*=1.11
@test worst_relative_error(out_temp,out_pulse)>1e-1
out_temp = deepcopy(out)
out_temp.T[1]*=1.11
@test worst_relative_error(out_temp,out_pulse)>1e-1
out_temp = deepcopy(out)
out_temp.T[1]*=1.11
@test worst_relative_error(out_temp,out_pulse)>1e-1

# compare to a pulse output with bigger timesteps, adapative solving should make this work
out_for_resample = rk8(12000,1e-7, tes, 1000,0);
out_pulse_ts = pulse(1200,1e-6, tes, 1000,0);
out_ts = ModelTES.TESRecord(out_for_resample.T[1:10:end], out_for_resample.I[1:10:end], out_for_resample.R[1:10:end], 1e-6)
@test worst_relative_error(out_ts,out_pulse_ts)<1e-5
out_pulse_ts2 = pulse(120,1e-5, tes, 1000,0);
out_ts2 = ModelTES.TESRecord(out_for_resample.T[1:100:end], out_for_resample.I[1:100:end], out_for_resample.R[1:100:end], 1e-5)
@test worst_relative_error(out_ts2,out_pulse_ts2)<1e-5

# Integrate a pulse with 12000 samples, 1e-7 second spacing, 1000 eV energy, 2000 presamples from the higher biased version of the same tes
out2 = pulse(12000,1e-7, tes2, 1000, 2000);

# Calculate a stochastic noise 1000 eV pulse with 12000 samples and 2000 presmples
outstochastic = stochastic(12000,1e-7, tes, 1000,2000);

# many pulses in one trace
outmany = ModelTES.pulses(12000,1e-7, tes2, [1000,1000,2000,3000,1000,500,2000], collect(1:7)*2e-4);
@test length(times(outmany))==length(outmany.I)
#make the pulses arrive halfway between time points
outmany2 = ModelTES.pulses(12000,1e-7, tes2, [1000,1000,2000,3000,1000,500,2000], 0.5e-7+collect(1:7)*2e-4);
@test times(outmany)==times(outmany2)
# compare the difference between when the pulses arrive half way between time points, and 1 time point apart
# the integrated difference should be about a factor of two apart
a=sum(abs(outmany.I[2:end-1]-outmany2.I[2:end-1])) # pulses off by half a sample
b=sum(abs(outmany.I[2:end]-outmany.I[1:end-1])) # pulses off by one sample
@test isapprox(a,b/2,rtol=1e-2,atol=1e-5)


# Make the other TESs in tes_models.jl
ModelTES.lowEpix()
ModelTES.highEpix()
ModelTES.pholmes(50e-9)
tesL = ModelTES.LCLSII(100e-9)

# Compute the noise spectrum etc.
nmodel = NoiseModel(tesL, 1e-6)
freq = linspace(0, 5e5, 50)
psd = model_psd(nmodel, freq)
covar = model_covariance(nmodel, 20)


include("ihtes.jl")
