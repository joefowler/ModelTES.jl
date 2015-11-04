using Polynomials

"Returns (noise, A,B,C,D) where `noise` is a noise power spectral
density in A^2 per Hz, at each input frequency (Hz) given by `freq`.
A,B,C,D are the components of the noise, such that A+B+C+D=noise.
In order, they are TES Johnson noise; load resistor Johnson noise;
thermal fluctuation noise; and amplifier noise.

You can input the amplifier current noise `SI_amp`, in A^2/Hz, or you
can take the default value.
"
function noise(tes::IrwinHiltonTES, freq::Vector{Float64}, SI_amp=5e-22)
    const kB = 1.38e-23 # Boltzmann
    const F  = 1 # this term goes from 0 to one and depends on wether the thermal conductivity is ballaistic or diffusive, hardcoded as 1 for now
    SP_TFN = 4*kB*tes.T0^2*tes.G0*F
    SV_TES = 4*kB*tes.T0*tes.R0*(1+2*tes.beta) # TES voltage noise
    SV_L   = 4*kB*tes.T0*tes.Rl  # Load voltage noise
    # SI_amp = 5e-22 # Amps^2/Hz

    omega = 2*pi*freq  # Radians / sec
    sIomeg = (1-tes.tauplus/tes.taucc)*(1-tes.tauminus/tes.taucc)./((1+1im*omega*tes.tauplus).*(1+1im*omega*tes.tauminus)) /(tes.I0*tes.R0*(2+tes.beta))
    sIomeg2 = abs2(sIomeg)

    Inoise_TFN = sIomeg2*SP_TFN
    Inoise_amp = SI_amp
    Inoise_TES = SV_TES*tes.I0^2/tes.loopgain^2 * (1+(tes.tauthermal*omega).^2) .* sIomeg2
    Inoise_load = SV_L*tes.I0^2*(tes.loopgain-1)^2/tes.loopgain^2 * (1+(tes.taucc*omega).^2) .* sIomeg2
    Inoise = Inoise_TFN+Inoise_amp+Inoise_TES+Inoise_load
    Inoise, Inoise_TES, Inoise_load, Inoise_TFN, Inoise_amp+zeros(Float64, length(Inoise))
end


"Returns noise of an IrwinHiltonTES modeled as an ARMA(2,2) process as
(theta,phi,sigma).

`theta` and `phi` are the coefficients of the MA(2) and AR(2) polynomials
to be used for approximating the noise as an ARMA(2,2) process, and `sigma`
is the rms current (Amps)."
function ARMAmodel(tes::IrwinHiltonTES, sampleTime::Float64, SI_amp=5e-22)
    const kB = 1.38e-23 # Boltzmann
    const F  = 1 # this term goes from 0 to one and depends on wether the thermal conductivity is ballaistic or diffusive, hardcoded as 1 for now
    SP_TFN = 4*kB*tes.T0^2*tes.G0*F
    SV_TES = 4*kB*tes.T0*tes.R0*(1+2*tes.beta) # TES voltage noise
    SV_L   = 4*kB*tes.T0*tes.Rl  # Load voltage noise

    A = SI_amp
    B = abs2((1-tes.tauplus/tes.taucc)*(1-tes.tauminus/tes.taucc) /(tes.I0*tes.R0*(2+tes.beta)))
    C = SP_TFN
    D = SV_TES*tes.I0^2/tes.loopgain^2
    E = SV_L*tes.I0^2*(tes.loopgain-1)^2/tes.loopgain^2

    # (u,v) are the solutions to the quadratic equation ax^2+bx+c=0
    a = A * (tes.tauplus*tes.tauminus)^2
    b = -(A*(tes.tauplus^2+tes.tauminus^2)+B*(D*tes.tauthermal^2+E*tes.taucc^2))
    c = A+B*(C+D+E)
    discr = sqrt(b^2-4a*c)
    u2 = (-b-discr)/(2a)
    v2 = (-b+discr)/(2a)

    u,v = sqrt(u2), sqrt(v2)
    K = sqrt(a)
    theta = [(sampleTime*u*.5+1)*(sampleTime*v*0.5+1),
            sampleTime^2*u*v*.5-2,
            (sampleTime*u*.5-1)*(sampleTime*v*.5-1)]
    phi = [(2tes.tauplus/sampleTime+1)*(2tes.tauminus/sampleTime+1),
            2-8tes.tauplus*tes.tauminus/sampleTime^2,
            (1-2tes.tauplus/sampleTime)*(1-2tes.tauminus/sampleTime)]

    # Rescale so that theta and phi have unit coefficients for the z^0 terms.
    sigma = theta[1]/phi[1]*4K/(sampleTime^2)
    theta ./= theta[1]
    phi ./= phi[1]

    # Note that  2sigma/sqrt(2pi*sampleTime) is the "sigma" needed in ARMA noise generation.
    theta, phi, sigma*2/sqrt(2pi*sampleTime)
end
ARMAmodel(tes::BiasedTES, sampleTime::Float64, SI_amp::Float64=5e-22) =
    ARMAmodel(IrwinHiltonTES(tes), sampleTime, SI_amp)

"Compute the noise power spectrum for an ARMA model with MA coefficients `theta`,
AR coefficients `phi`, and overall noise scale of `sigma.`  Use the sample time `T`
(seconds) and compute PSD at (natural, not angular) frequencies `freq` (Hz).
Note that `p1` and `p2` should be roughly the same in the following example:

SI_amp = 8e-22
freq = logspace(1,6,50)
T = 0.5/freq[end]
(theta,phi,sigma) = ARMAmodel(tes, T, SI_amp)
p1 = ARMApowerspectrum(theta, phi, sigma, freq, T)
p2 = noise(tes, freq, SI_amp)
"
function ARMApowerspectrum(theta::Vector, phi::Vector, sigma::Float64, freq::Vector, T::Float64)
    # Normalize the coefficients
    sigma *= theta[1]/phi[1]
    theta = theta / theta[1]
    phi = phi / phi[1]

    omega_digital = 2*pi*freq*T # Proportional to freq, and =pi at the Nyquist frequency
    z = exp(-1im*omega_digital)
    Inoise = sigma^2*abs2((1+theta[2]*z+theta[3]*(z .^2)) ./ (1+phi[2]*z+phi[3]*(z .^2)))
end



"Compute a covariance vector for the ARMA(p,q) process represented
by the polynomials `theta` for the MA(q) and `phi` for the AR(p) part
of the process. `sigma` is the digital version of the noise scale, in
units of output per sqrt(Hz). `sampleTime` is the time between samples
in seconds. `Nsamp` values will be computed and returned"

function ARMAcovariance(theta::Vector, phi::Vector, sigma::Float64, sampleTime::Float64, Nsamp::Integer)
    k, xi, covar0 = analyzeARMA22(theta, phi, sigma)
    t = collect(0:-1:1-Nsamp)
    covar = k[1]*(xi[1].^t) + k[2]*(xi[2].^t)
    covar[1] = covar0
    covar
end


"Find numbers needed to compute the covariance vector for the ARMA(p,q) process
given by the polynomials `theta` for the MA(q) and `phi` for the AR(p) part
of the process. `sigma` is the digital version of the noise scale, in
units of output per sqrt(Hz).
Returns (`b`,`xi`,`c0`), where `c0` is the covariance at lag 0, and for any
lag j>0, the covariance is the sum over i={1,2} of b[i]*xi[i]^j.
"

function analyzeARMA22(theta::Vector{Float64}, phi::Vector{Float64}, sigma::Number)
    # Only works on ARMA(2,2) processes. We know how to generalize, but why? The
    # Irwin-Hilton model guarantees ARMA(2,2) noise model.
    p,q = length(phi)-1, length(theta)-1
    @assert p == 2
    @assert q == 2
    psi = [1, theta[2]-phi[2], theta[3]-phi[3]-phi[2]*theta[2]+phi[2]^2]
    rhs = sigma^2 * [psi[1]+theta[2]*psi[2]+theta[3]*psi[3],
                theta[2]*psi[1]+theta[3]*psi[2], theta[3]*psi[1]]
    M = [1 phi[2] phi[3];
        phi[2] 1+phi[3] 0;
        phi[3] phi[2] 1]
    gamma = M \ rhs

    pm = sqrt(phi[2]^2-4*phi[1]*phi[3])
    xi = ([pm, -pm]-phi[2])*(0.5/phi[3]) # roots of phi polynomial a1+a2z+a3z^2
    M2 = [xi' .^ -1;  xi' .^ -2]
    b = M2 \ gamma[2:3]
    # How to compute the full covariance:
    # covar=zeros(Float64, 1000)
    # covar[1:3] = gamma[1:3]
    # xiterms = 1. ./ xi
    # for i = 2:length(covar)
    #     covar[i] = dot(xiterms, b)
    #     xiterms ./= xi
    # end
    b, xi, gamma[1]
end


"""Generate N samples of ARMA noise, where
`theta` is the vector of coefficients for the MA polynmial of order p
`phi` is the vector of coefficients for the AR polynmial of order q
`sigma` is the the sqrt of the output variance
"""

function generateARMAnoise(theta::Vector{Float64}, phi::Vector{Float64}, sigma::Number, N::Int)
    p,q = length(phi)-1, length(theta)-1
    # Rescale the rational function to have 0-order terms of unit size.
    if theta[1] != 1
        sigma *= theta[1]
        theta = theta ./ theta[1]
    end
    if phi[1] != 1
        sigma /= phi[1]
        phi = phi ./ phi[1]
    end

    # Need to discard data for MA-related transients to fade out
    # Let the discarded data be at least 20 e-folding times for the
    # slowest root of theta(z), OR 10k samples, whichever is longer.
    r = Polynomials.roots(Poly(theta, :x))
    fadetime = 20./log(r)
    const Ndiscard=max(10000,Int(floor(maximum(fadetime))))
    const Ncompute = N+Ndiscard

    # MA part of the process
    eps = randn(Ncompute) * sigma
    out = copy(eps)
    for i = q+1:Ncompute
        for j=1:q
            out[i] += eps[i-j]*theta[j+1]
        end
    end

    # AR part of the process
    for i = 1+p:Ncompute
        for j = 1:p
            out[i] -= out[i-j]*phi[j+1]
        end
    end
    out[Ndiscard+1:end]
end
