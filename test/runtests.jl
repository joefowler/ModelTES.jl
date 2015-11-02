using ModelTES
using Base.Test

# write your own tests here
tes = BiasedTES(Holmes48nH, Holmes48nH.Rn*0.2)
tes2 = BiasedTES(Holmes48nH, Holmes48nH.Rn*0.4)
out = rk8(3000,1e-6, tes, 1000)
out2 = rk8(3000,1e-6, tes2, 1000)
linearparams = getlinearparams(tes)
lintes = IrwinHiltonTES(linearparams...)
f = logspace(0,6,100)
n,n1,n2,n3 = noise(lintes, f)
# figure()
# plot(times(out), out.I)
# plot(times(out2), out2.I)
# figure()
# loglog(f,n)
# loglog(f,n1,label="ites")
# loglog(f,n2,label="il")
# loglog(f,n3,label="itfn")
# legend()
