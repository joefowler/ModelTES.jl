using ModelTES, PyPlot, Formatting

bt = ModelTES.lowEpix()

"biaspulse(bt, Vratio, dt, npre, nheat, noff)
Return a pulse record with `dt*npre` seconds of pretrigger, followed by
`dt*nheat` seconds of the bias voltage being multiplied by `Vratio` (heating if `Vratio`>1)
followed by `dt*noff` seconds of relaxation (bias voltage returned to normal value = bt.V)."
function biaspulse(bt, Vratio, dt, npre, nheat, noff)
  bt2 = ModelTES.BiasedTES(bt.p, bt.I0, bt.T0, bt.V*Vratio)
  #function pulse(nsample::Int, dt::Float64, bt::BiasedTES, E::Number, npresamples::Int=0
  pre = ModelTES.TESRecord(fill(bt.T0,npre), fill(bt.I0, npre), fill(ModelTES.getR0(bt),npre), dt)
  heat_on = pulse(nheat+1,dt, bt2, 0, 0)
  bt3 = ModelTES.BiasedTES(bt.p, heat_on.I[end], heat_on.T[end], bt.V)
  heat_off = pulse(noff+1, dt, bt3,0, 0)

  @assert pre.I[end]==heat_on.I[1]
  @assert pre.T[end]==heat_on.T[1]
  @assert heat_on.I[end]==heat_off.I[1]
  @assert heat_on.T[end]==heat_off.T[1]

  jointI = cat(1,pre.I, heat_on.I[2:end], heat_off.I[2:end])
  jointT = cat(1,pre.T, heat_on.T[2:end], heat_off.T[2:end])

  joint = ModelTES.TESRecord(jointT, jointI, [ModelTES.R(jointI[i],jointT[i],bt.p) for i=1:length(jointT)],dt)
  @assert length(joint.I)==npre+nheat+noff
  return joint
end


Vratio = 1.5
dt = 1e-7
npre = 10000
nheat = 50000
dtheat = dt*nheat
noff = 50000
joint = biaspulse(bt, Vratio, dt, npre, nheat, noff)
Piv = [joint.I[i]^2*joint.R[i] for i=1:length(joint.I)]
Plegs = bt.p.k*(joint.T.^bt.p.n-bt.p.Tbath^bt.p.n)
V = cat(1,fill(bt.V,npre), fill(bt.V*Vratio, nheat), fill(bt.V, noff))
Rshunt = bt.p.Rl-bt.p.Rp
Ibias = V/Rshunt
# Ereffactor = 0.9715
Eref0 = (maximum(joint.T)-joint.T[1])*bt.p.C/ModelTES.J_per_eV
pulseref0 = pulse(length(joint.I), dt, bt, Eref0, 0)
j = indmin(joint.I)-indmin(pulseref0.I)
Ereffactor = (maximum(joint.T)-bt.T0)/(maximum(pulseref0.T)-bt.T0)
Eref = Eref0*Ereffactor
pulseref = pulse(length(joint.I), dt, bt, Eref , j)

figure(figsize=(16,10))
subplot(221)
plot(times(joint),joint.I,label="bias pulse",lw=2)
plot(times(pulseref),pulseref.I,label="$(format(Eref, precision=2)) eV",lw=2)
xlabel("time (s)")
ylabel("TES current (A)")
legend(loc="best")

subplot(222)
plot(times(joint),joint.T,label="bias pulse",lw=2)
plot(times(pulseref),pulseref.T,label="$(format(Eref, precision=2))  eV",lw=2)
xlabel("time (s)")
ylabel("TES temp (K)")
legend(loc="best")

subplot(223)
plot(times(joint), Piv,lw=2, label="IV")
plot(times(joint), Plegs,lw=2, label="legs")
xlabel("time (s)")
ylabel("TES power (W)")
legend(loc="best")

subplot(224)
plot(times(joint), Ibias, lw=2)
xlabel("time (s)")
ylabel("Room temperatue bias current (A)")


#calculate pulseheight change vs Tbath for bias pulse and normal pulse
dT = 0.001
Tbs = linspace(bt.p.Tbath-dT, bt.p.Tbath+dT,100)
function newbt(Tbath, bt)
  bt0 = deepcopy(bt)
  bt0.p.Tbath = Tbath
  r=pulse(10,1e-1, bt0, 0)
  return ModelTES.BiasedTES(bt0.p, r.I[end], r.T[end], bt0.V)
end
bts = [newbt(Tbath, bt) for Tbath in Tbs]
@assert all([bt.V for bt in bts].==bt.V)

bias_pulses = [biaspulse(bt, Vratio, dt, npre, nheat, noff) for bt in bts]
normal_pulses = [pulse(length(joint.I), dt, bt, Eref, j) for bt in bts]
ph_bias_pulses = [minimum(record.I)-record.I[1] for record in bias_pulses]
ph_normal_pulses = [minimum(record.I)-record.I[1] for record in normal_pulses]

ph_bias_ref = minimum(joint.I)-joint.I[1]
ph_normal_ref = minimum(pulseref.I)-pulseref.I[1]
figure()
plot(Tbs, ph_bias_pulses./ph_bias_ref,label="bias pulse", lw=2)
plot(Tbs, ph_normal_pulses./ph_normal_ref,label="$(format(Eref, precision=2))  eV",lw=2)
xlabel("Tbath (K)")
ylabel("fractional peakheight change")
legend(loc="best")
