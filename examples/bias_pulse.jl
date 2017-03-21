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
P = [joint.I[i]^2*joint.R[i] for i=1:length(joint.I)]
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
plot(times(joint),joint.I,label="heater",lw=2)
plot(times(pulseref),pulseref.I,label="$(format(Eref, precision=2)) eV",lw=2)
xlabel("time (s)")
ylabel("TES current (A)")
legend(loc="best")

subplot(222)
plot(times(joint),joint.T,label="heater",lw=2)
plot(times(pulseref),pulseref.T,label="$(format(Eref, precision=2))  eV",lw=2)
xlabel("time (s)")
ylabel("TES temp (K)")
legend(loc="best")

subplot(223)
plot(times(joint), P,lw=2)
xlabel("time (s)")
ylabel("TES IV power (W)")

subplot(224)
plot(times(joint), Ibias, lw=2)
xlabel("time (s)")
ylabel("Room temperatue bias current (A)")
