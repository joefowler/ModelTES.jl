using ModelTES
using Base.Test

Tc=0.107; Tbath = 0.065
R0 = 1.55e-3; Rl = 0.35e-3; Rn = 10.3e-3; Rpara=0.0
n=3.3; k = 3.50e-9; C = 0.72e-12
alpha = 92.0;β = 1.4; L=300e-9
model = ModelTES.ShankRIT(alpha, β, n, Tc, Tbath, k, R0, Rn);
tes_param = TESParams(n,Tc,Tbath,k,C,L,Rl,Rpara,Rn,model)
bt = BiasedTES(tes_param, R0)
ihtes = IrwinHiltonTES(bt)

@test Tbath == ihtes.Tbath
@test isapprox(R0, ihtes.R0; rtol=1e-5)
@test Rl == ihtes.Rl
@test C == ihtes.C0
@test isapprox(alpha, ihtes.alpha; rtol=1e-5)
@test isapprox(β, ihtes.beta; rtol=1e-5)
@test L == ihtes.L

f = logspace(0,6,100);
n,n1,n2,n3,n4 = ModelTES.noisePSD(ihtes, f);
z = ModelTES.Z(ihtes,f)
zcircuit = ModelTES.Zcircuit(ihtes,f)

zr, zi = real.(z), imag.(z)
