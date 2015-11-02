using ModelTES
using Base.Test

# createa a Biased TES from the 48 nanohentry Holmes paramters with 0.2*Rn resistance
tes = BiasedTES(Holmes48nH, Holmes48nH.Rn*0.2)
# createa a Biased TES from the 48 nanohentry Holmes paramters with 0.4*Rn resistance
tes2 = BiasedTES(Holmes48nH, Holmes48nH.Rn*0.4)
# integrate a pulse with 12000 samples, 1e-7 second spacing, 1000 eV energy, 2000 presamples
out = rk8(12000,1e-7, tes, 1000, 2000)
# integrate a pulse with 12000 samples, 1e-7 second spacing, 1000 eV energy, 2000 presamples from the higher biased version of the same tes
out2 = rk8(12000,1e-7, tes2, 1000, 2000)
# get all the linear parameters for the irwin hilton model
linearparams = getlinearparams(tes)
# store them in an IrwinHiltonTES type
lintes = IrwinHiltonTES(linearparams...)
# caulculate the noise and the 3 component terms in the IrwinHilton model
f = logspace(0,6,100)
n,n1,n2,n3 = noise(lintes, f)


#calculate a stocastic noise 1000 eV pulse with 12000 samples and 2000 presmples
outstochastic = stochastic(12000,1e-7, tes, 1000,2000)

# using PyPlot
# figure()
# same parameters, one noiseless, one with stocastic noise
# plot(times(out),out.I,".")
# plot(times(outstochastic),outstochastic.I,".")
# figure()
# title("pulses with different bias points")
# plot(times(out), out.I)
# plot(times(out2), out2.I)
# figure()
# loglog(f,n)
# loglog(f,n1,label="ites")
# loglog(f,n2,label="il")
# loglog(f,n3,label="itfn")
# legend()
