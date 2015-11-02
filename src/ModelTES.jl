module ModelTES
export transitionwidth, BiasedTES, TESParams, getlinearparams, noise, noiseterms,
TESRecord, times, rk8, Holmes48nH, IrwinHiltonTES, stochastic
using Roots, ForwardDiff
include("rk8.jl")


abstract AbstractRIT
# following Irwin-Hilton figure 3
type TESParams{RITType<:AbstractRIT}
    n       ::Float64   # thermal conductance exponent (unitless)
    Tc      ::Float64   # TES critical temperature (K)
    Tbath   ::Float64   # bath temperature (K)

    k       ::Float64   # Pbath pre-factor (W/K^n)
    C       ::Float64   # heat capacity of TES (J/K)

    L       ::Float64   # inductance of SQUID (H)
    Rl      ::Float64   # Thevenin-equivalent resistance Rload = (Rshunt+Rparasitic)(ohms)
    Rp      ::Float64   # Rparastic, Rshunt = Rl-Rparasitic (ohms)
    Rn      ::Float64   # normal resistance of TES (ohms)

    RIT     ::RITType   # RIT surface

end
transitionwidth(p::TESParams)=transitionwidth(p.RIT)


type ShankRIT <: AbstractRIT
    Tw::Float64 # transition width (K)
    A::Float64  # current dependence for R(T,I) (A/K^(3/2))
end
transitionwidth(RIT::ShankRIT)=RIT.Tw
"Calculate `Tw` and `A` for a ShankRIT to have the given `alpha` and `beta` paramters
when biased at resistane `R0`."
function ShankRIT(alpha, beta, n, Tc, Tb, k, R0, Rn)
    T0 = Tc /(1 + 3*beta/(2*alpha) - 2*(Rn-R0)/(Rn*alpha)*atanh(2*R0/Rn-1))
    I0 = sqrt(k*(T0^n-Tbath^n)/R0)
    Tw = T0*(Rn-R0)/(Rn*log(2)*alpha)
    A = I0*(2*alpha/(3*T0*beta))^(3/2)
    ShankRIT(Tw,A)
end

type BiasedTES{T}
    p::TESParams{T}
    I0::Float64 # intial current for diff eq (A)
    T0::Float64 # initial temperature for diff equations (K)
    V ::Float64 # thevinen equivalent voltage V = I0*(p.Rl+p.R0), R0=quiescent resistance
end
"Used to find the value of T0 that makes R==targetr (or at least come close)."
function getR0error(T0,targetr, p::TESParams)
   I0 = sqrt(p.k*(T0^p.n-p.Tbath^p.n)/targetr)
   R(I0,T0,p)-targetr
end
"Find the initial conditions (I0, T0, V) that cause `p` to have resistance `targetr`."
function initialconditions(p::TESParams, targetr)
   T00 = fzero((t)->getR0error(t, targetr,p),p.Tc-10*transitionwidth(p), p.Tc+10*transitionwidth(p))
   I00 = sqrt(p.k*(T00^p.n-p.Tbath^p.n)/targetr)
   R00 = R(I00,T00,p)
   V00 = I00*(p.Rl+R00)
   # now evolve these conditions through integration to really lock themin
   out = rk8(1000,1e-5, BiasedTES(p, I00, T00, V00), 0)
   T0 = out.T[end]
   I0 = out.I[end]
   R0 = R(I0,T0,p)
   V = I0*(p.Rl+R0)
   I0,T0,V
end
"Created a biased tes with quiescent state resistance R0"
function BiasedTES{T<:ShankRIT}(p::TESParams{T}, R0::Float64)
   @assert 0<R0<p.Rn
   I0, T0, V = initialconditions(p,R0)
   BiasedTES(p,I0,T0,V)
end
"Calculate `R0` the quiescent resistance of `bt`."
getR0(bt::BiasedTES) = R(bt.I0, bt.T0, bt.p)
"Calculate thermal conductivity from TES heat capacit to Tbath in small signal limit."
function getG0(bt::BiasedTES)
   p=bt.p
   p.k*p.n*bt.T0^(p.n-1)
end
"Calculate paramaters in Irwin-Hilton table 1."
function getlinearparams(bt::BiasedTES)
   f(x) = R(x[1], x[2], bt.p)
   p = bt.p
   R0 = getR0(bt)
   G0 = getG0(bt)
   g = ForwardDiff.gradient(f)
   drdi,drdt = g([bt.I0, bt.T0])
   alpha = drdt*bt.T0/R0
   beta = drdi*bt.I0/R0
   PJ = bt.I0^2*R0
   loopgain = PJ*alpha/G0/bt.T0
   tauthermal = p.C/G0
   taucc = tauthermal/(1-loopgain) # constant current time constant
   r = p.Rl/R0
   taueff = (1+beta+r)/(1+beta+r+(1-r)*loopgain) # zero inductance effecticve thermal time constant
   tauelectrical = p.L/(p.Rl+R0*(1+beta))
   invtau = 1/(2*tauelectrical)+1/(2*taucc)
   a = (1/tauelectrical-1/taucc)^2
   b = -4*(R0/p.L)*loopgain*(2+beta)/tauthermal
   invtaupm = 0.5*sqrt(a+b)
   tauplus = 1/(invtau+invtaupm)
   tauminus = 1/(invtau-invtaupm)
   c = loopgain*(3+beta-r)+(1+beta+r)
   d = loopgain*(2+beta)*(loopgain*(1-r)+(1+beta+r))
   f = 2*sqrt(d)*R0*tauthermal/(loopgain-1)^2
   lcritplus = c+f
   lcritminus = c-f
   bt.I0, bt.T0, bt.V, p.Rl, p.Tbath, p.Tbath, p.L, R0, G0, p.C, alpha, beta, loopgain, tauthermal, taucc, taueff, tauelectrical, tauplus, tauminus, lcritplus, lcritminus
end

"Paramters from Irwin-Hilton table one for modeling a linear TES."
type IrwinHiltonTES
   I0::Float64
   T0::Float64
   V::Float64
   Rl::Float64
   Tl::Float64 # temperature of the load resistor, usually modeled as =Tbath, but really should come from modeling ep coupling in the load resistor
   Tbath::Float64
   L::Float64
   R0::Float64
   G0::Float64
   C0::Float64
   alpha::Float64
   beta::Float64
   loopgain::Float64
   tauthermal::Float64
   taucc::Float64
   taueff::Float64
   tauelectrical::Float64
   tauplus::Float64
   tauminus::Float64
   lcritplus::Float64
   lcritminus::Float64
end

isoverdamped(tes::IrwinHiltonTES) = isreal(tes.tauplus) && isreal(tes.tauminus) && tes.tauplus<tes.tauminus
isunderdamped(tes::IrwinHiltonTES) = !isoverdamped(tes) && tes.tauplus!=tes.tauminus
function powertocurrent(I0, T0, V, Rl, Tl, Tbath, L, R0, G0, C, alpha, beta, loopgain, tauthermal, taucc, taueff, tauelectrical, tauplus, tauminus, lcritplus, lcritminus, omega)
   a = -1/I0/R0
   b = L/tauelectrical/R0/loopgain + (1-Rl/R0)
   c = im*L*tauthermal/R0/loopgain*(1/taucc+1/tauelectrical)
   d = -tauthermal*L/loopgain/R0
   -a.*(b.+omega.*c.+omega.^2.*d).^-1
end
function noiseterms(I0, T0, V, Rl, Tl, Tbath, L, R0, G0, C, alpha, beta, loopgain, tauthermal, taucc, taueff, tauelectrical, tauplus, tauminus, lcritplus, lcritminus, omega)
   xi = 1+2*beta #quadratic approximation
   F  = 1 # this term goes from 0 to one and depends on wether the thermal conductivity is ballaistic or diffusive, hardcoded as 1 for now
   si = abs2(powertocurrent(I0, T0, V, Rl, Tl, Tbath, L, R0, G0, C, alpha, beta, loopgain, tauthermal, taucc, taueff, tauelectrical, tauplus, tauminus, lcritplus, lcritminus, omega))
   SItes = 4*kb*T0*I0^2*R0*xi/loopgain^2*(1+omega.^2*tauthermal).*si
   SIl   = 4*kb*Tl*I0^2*Rl*(loopgain-1)^2/loopgain^2*(1+omega.^2*taucc^2).*si
   SItfn = 4*kb*T0^2*G0*F*si
   SItes, SIl, SItfn
end
function noise(tes::IrwinHiltonTES, f)
   omega = 2*pi*f
   args = Any[getfield(tes, fieldname) for fieldname in fieldnames(tes)]
   push!(args,omega)
   SItes, SIl, SItfn = noiseterms(args...)
   SItes+SIl+SItfn, SItes, SIl, SItfn
end


type TESRecord
    T::Vector{Float64}  # temperature (K)
    I::Vector{Float64}  # current (A)
    dt::Float64   # seconds between samples (seconds)
end
times(r::TESRecord) = range(0,r.dt,length(r.I))
# conversion J/eV
const J_per_eV = 1.602177e-19 #unitless
const kb = 1.38064852e-23 #k boltzmann (J/K)

# I want to re-write this code so that
# A and Tw are the TES params, alpha and beta are derived
# and R0 and be changed to change the bias conditions
# then I0 and T0 and derived from those 3 parametrs

"find bias parameters T0, I0, V for to make the TES have resistane R0 in steady state"
function findbiasparams(tes::TESParams{ShankRIT},R0)

    T0 = Tc /(1 + 3*beta/(2*alpha) - 2*(Rn-R0)/(Rn*alpha)*atanh(2*R0/Rn-1))

    # V = I0*(p.Rl+p.R0)
end

"TES resistance model (Shank et al. 2014)"
R(I, T, RIT::ShankRIT, Tc, Rn) = Rn/2*(1+tanh((T-Tc+(max(I,0.0)/RIT.A)^(2/3))/(2*log(2)*RIT.Tw)))
R(I,T, p::TESParams) = R(I,T, p.RIT, p.Tc, p.Rn)
"thermal TES equation
C*dT/dt = I^2*R - k*(T^n-Tbath^n)"
function dT(I, T, k, n, Tbath, C, R)
    Q=-k*(T^n-Tbath^n)+I^2*R
    Q/C
end

"electrical TES equation
L*di/dt = (IBias-I)*Rs-I*Rs-I*Rtes"
function dI(I, T, V, Rl, L, R)
    (V-I*(Rl+R))/L
end

"Modify `dy` to contain "
function stochasticdy!(dy,y,dt,bt::BiasedTES)
   p=bt.p
   T, I = y # tes tempertature, tes current
   r = R(I,T,p) #res resistance
   # generate random variables for use in stocastic integration
   bath_thermal_noise, tes_thermal_noise, rl_johnson_noise, tes_johnson_noise = randn(4)

   johnsontes = tes_johnson_noise*2*sqrt(kb*T*r/dt)
   johnsonrl = rl_johnson_noise*2*sqrt(kb*p.Tbath*p.Rl/dt)

   # k*T^n*dt/kb*T is the average number of phonons going from the TES to the bath
   # assuming all phonons of energy kb*t
   # sqrt(k*T^n*dt/kb*T) is the std deviation of the number of phonons
   # the same reasoning applies in the other direction
   dqdt = -p.k*(T^p.n-p.Tbath^p.n)+I^2*r +
     tes_thermal_noise*sqrt(kb*p.k/dt)*T^(p.n/2 + 0.5) + # variation in number of phonon TES->bath
     bath_thermal_noise*sqrt(kb*p.k/dt)*p.Tbath^(p.n/2+0.5) + # variation in number of phonon bath->TES
     -I*johnsontes # work on the noise voltage source for the TES add/remove energy from the TES

   # the resistors are modeled as noiseless resistors in parallel with
   # voltage sources having white noise
   # the load resistor is assumed to be at Tbath
   ldidt = bt.V-I*(p.Rl+r) + johnsontes + johnsonrl

   if isnan(dqdt)
      @show dy, y, dt, bt
      error()
   end


   dy[1] = dt*dqdt/p.C
   dy[2] = dt*ldidt/p.L
end

"stochastic(nsample::Int, dt::Float64, bt::BiasedTES, E::Number, npresample::Int=0, nrandomizingsample::Int=6000)
Integrate a pulse record with with stocastic noise `E` eV with `nsample` total samples, `npresample` presamples, and use
`nrandomizingsample` samples thrown away to ensure randomized starting `T0` and `I0`. Use `E=0` for noise records."
function stochastic(nsample::Int, dt::Float64, bt, E::Number, npresample::Int=0, nrandomizingsample::Int=6000)
   p = bt.p
   TI = Array(Float64, 2, nsample)
   y = [bt.T0, bt.I0]
   dy = Array(Float64,2)
   # manually implement euler's method for integration
   for i = 1:nrandomizingsample # integrate for a while to get well randomized T0 and I0
       stochasticdy!(dy,y,dt,bt)
       y+=dy
   end
   for i = 1:nsample
      if i == npresample+1
         y[1]+=E*J_per_eV/p.C
      end # increase temperature to start pulse
      stochasticdy!(dy,y,dt,bt)
      y+=dy
      TI[:,i]=y
   end
   TESRecord(reshape(TI[1,:],(nsample,)),reshape(TI[2,:],(nsample,)),dt)
end

"Calling a BiasedTES gives the dI and dT terms for integration."
function Base.call{S<:Float64}(bt::BiasedTES, t::Float64, y::AbstractVector{S}, dy::AbstractVector{S})
    T,I = y[1],y[2]
    p = bt.p
    r = R(I,T,p)
    dT(I, T, p.k, p.n, p.Tbath, p.C, r)
    dI(I,T, bt.V, p.Rl, p.L, r)
    dy[:] = [dT(I, T, p.k, p.n, p.Tbath, p.C, r),
             dI(I,T, bt.V, p.Rl, p.L, r)]
end

function rk8(nsample::Int, dt::Float64, bt::BiasedTES, E::Vector{Float64}, npresamples::Int=0)
    out = Vector{TESRecord}(length(E))
    for i in 1:length(E)
        out[i] = rk8(nstep, dt, bt, E[i], npresamples)
    end
    out
end

function rk8(nsample::Int, dt::Float64, bt::BiasedTES, E::Number, npresamples::Int=0)
    # Pair of differential equations y' = f(t,y), where y=[T,I]
    p = bt.p
    # Integrate pair of ODEs for all energies EE
    TI = Array(Float64, 2, nsample)
    y = Array(Float64, 2); ys = similar(y); work = Array(Float64, 14)
    TI[1,:] = bt.T0+E*J_per_eV/p.C # set T0
    TI[2,:] = bt.I0 # set I0
    y[:] = TI[:,1]
    for i = npresamples+1:nsample
        rk8!(bt, 0.0, dt, y, ys, work)
        TI[:, i] = ys
        y[:] = ys
    end
    TESRecord(reshape(TI[1,:],(nsample,)), reshape(TI[2,:], (nsample,)), dt)
end

Holmes48nH = TESParams(3.25, 0.1, 0.07, 2.33e-8, 0.5e-12, 48e-9, 0.3e-3, 0.0, 10e-3, ShankRIT(0.00056, 1.13))

end # module
