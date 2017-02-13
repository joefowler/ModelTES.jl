using ModelTES, Interpolations

type InterpolatedRIT{T} <: ModelTES.AbstractRIT
  interpolator::T
end

T = collect((100:0.1:106)/1000)
slopes = ones(length(T)-1)
slopes[1:2]=0
# slopes[4:5]=1.1
slopes[17:20]=1.3
R = cumsum(slopes.*diff(T))
push!(R,R[end])
R*=0.01/R[end]
a=interpolate((T,),R,Gridded(Linear()))



ModelTES.R(I, T, RIT::InterpolatedRIT, Tc, Rn) = RIT.interpolator[T,I]
ModelTES.transitionwidth(RIT::InterpolatedRIT) = std(T)
p = ModelTES.TESParams(
3.3, #n
mean(T), #Tc
0.07, #Tbath
3.5e-9, #k
10e-13, #C
1.5e-7, #L
0.00035, #Rl
0.0, #Rp
R[end], #Rn
InterpolatedRIT(a), #RIT
)

bt = BiasedTES(p, p.Rn*0.15)
ih = ModelTES.IrwinHiltonTES(bt)

function alpha(p::ModelTES.TESParams, T)
  f(x) = ModelTES.R(0, x[1], bt.p)
  r = ModelTES.R(0,T,bt.p)
  drdt = ForwardDiff.gradient(f, [T])
  alpha = drdt[1]*T/r
end

using PyPlot
Ts = linspace(bt.T0,maximum(T),100)
plot(Ts, [alpha(bt.p,t) for t in Ts])
plot(bt.T0, ih.alpha,"o")
xlabel("T (K)")
ylabel("α")
function pulsefb(energy, bt=bt)
  r=ModelTES.rk8(40000,1e-7,bt,energy,100)
  fb = -(r.I-r.I[1])
  fb
end
# function maxt(energy, bt=bt)
#   r=ModelTES.pulse(4000,1e-6,bt,energy,100)
#   maximum(r.T)
# end
energies = 1:100:10000
pulses = pulsefb.(energies)
# maxts = maxt.(energies)
tstarts = bt.T0+(energies./p.C*ModelTES.J_per_eV)
alphas = [alpha(bt.p,t) for t in tstarts]
av = mean(pulses)
fv = [dot(av,pulse) for pulse in pulses]

figure()
for j = 10:10:100
  plot(pulses[j])
end
""

using HDF5
h5 = h5open("out","w")
h5["energies"]=collect(energies)
h5["filt_value"]=fv
close(h5)
