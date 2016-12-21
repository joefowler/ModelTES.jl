using ModelTES, DifferentialEquations


tes = ModelTES.highEpix()

@time out = rk8(12000,1e-7, tes, 1000, 2000);

@time out2 = ModelTES.pulse(12000,1e-7, tes, 1000, 2000);


nsample=12000
dt = 1e-7
bt = tes
E=1000
npresamples = 2000
p = bt.p
# u = [T,I]
u0 = [bt.T0+E*ModelTES.J_per_eV/p.C, bt.I0]
function du(t,u)
    T,I = u[1],u[2]
    r = ModelTES.R(I,T,p)
    [ModelTES.dT(I, T, p.k, p.n, p.Tbath, p.C, r),
             ModelTES.dI(I,T, bt.V, p.Rl, p.L, r)]
end
prob = ODEProblem(tes, u0, (0.0,nsample*dt))
@time sol = DifferentialEquations.solve(prob,DifferentialEquations.Vern8(),dt=1e-9,abstol=1e-11,reltol=1e-11);
@time out =  sol(0:dt:(nsample-1)*dt);
@time sol2 = DifferentialEquations.solve(prob,DifferentialEquations.Vern8(),dt=1e-9,abstol=1e-11,reltol=1e-11, saveat=0:dt:(nsample-1)*dt, dense=false);


T=[o[1] for o in out];
I=[o[2] for o in out];
ModelTES.TESRecord(T,I, ModelTES.R(I,T,bt.p),dt)



function lorenz(t,u,du)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
lu0 = [1.0;0.0;0.0]
tspan = (0.0,1.0)
dt = 0.0001
times = first(tspan):dt:last(tspan)
method = DifferentialEquations.Tsit5()
lprob = ODEProblem(lorenz,lu0,tspan)
lt1 = @elapsed lsol = solve(lprob, method);
lt2 = @elapsed lsol(times);
lt3 = @elapsed lsol2 = solve(lprob, method, saveat=times, dense=false, save_timeseries=false);
lsol3 = solve(lprob, method, saveat=collect(times), dense=false, save_timeseries=false);
@show length(lsol2.t)
@show length(lsol3.t)

dt=1e-5
nsample=1000
prob = ODEProblem(tes, u0, (0.0,nsample*dt))
times2 = 0:dt:(nsample-1)*dt
t1 = @elapsed sol = solve(prob,method,dt=1e-9,reltol=1e-9);
t2=@elapsed out =  sol(0:dt:(nsample-1)*dt);
t3=@elapsed sol2 = solve(prob,method,dt=1e-9,reltol=1e-9, saveat=times2, dense=false);
t3=@elapsed sol2_abstol = solve(prob,method,dt=1e-9,abstol=1e-9,reltol=1e-9, saveat=times2, dense=false);

@printf("lorenz %s, solve %0.6f, sol interp %0.6f, solve saveat %0.6f\n", repr(typeof(method)), lt1, lt2, lt3);@printf("%s, solve %0.6f, sol interp %0.6f, solve saveat %0.6f", repr(typeof(method)), t1, t2, t3)


@time sol_saveat = solve(prob,method,dt=1e-9,reltol=1e-9, dense=false, saveat=times2, save_timeseries=false);

figure()
plot(sol2_abstol.t,sol2_abstol[:,2],label="abstol,reltol=1e-9")
plot(sol2.t,sol2[:,2],label="reltol=1e-9")
xlabel("time (s)")
ylabel("current (A)")
