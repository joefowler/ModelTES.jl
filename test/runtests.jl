using ModelTES
using Base.Test

using ModelTES
# Create a high-E TES design
stdTES = ModelTES.highEpix()
# Create a Biased TES from the 48 nanohentry Holmes paramters with 0.2*Rn resistance
tes = ModelTES.pholmes(48e-9, 0.20)
# Create a Biased TES from the 48 nanohentry Holmes paramters with 0.4*Rn resistance
tes2 = ModelTES.pholmes(48e-9, 0.40)
# Integrate a pulse with 12000 samples, 1e-7 second spacing, 1000 eV energy, 2000 presamples
out = rk8(12000,1e-7, tes, 1000, 2000);
# Integrate a pulse with 12000 samples, 1e-7 second spacing, 1000 eV energy, 2000 presamples from the higher biased version of the same tes
out2 = rk8(12000,1e-7, tes2, 1000, 2000);
# Get all the linear parameters for the irwin hilton model
lintes = IrwinHiltonTES(tes)
# Calculate the noise and the 4 components in the IrwinHilton model
f = logspace(0,6,100);
n,n1,n2,n3,n4 = noise(lintes, f);

# Calculate a stochastic noise 1000 eV pulse with 12000 samples and 2000 presmples
outstochastic = stochastic(12000,1e-7, tes, 1000,2000);
