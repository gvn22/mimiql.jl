using DifferentialEquations, LinearAlgebra
using Plots, PyPlot

# 2D GQL analogue: low mode & high mode
function lh!(du,u,p,t)

    x,y   = u
    α,β,γ = p

    du[1] = dx = α*x + β*(x*x + y*y) + γ
    du[2] = dy = α*y + β*(2.0*x*y)

end

# 2D GCE2 analogue: low mode & second cumulant
function lc!(du,u,p,t)

    x,z   = u
    α,β,γ = p

    du[1] = dx = α*x + β*(x*x + z) + γ
    du[2] = dz = 2.0*(α*z + β*(2.0*x*z))

end

# initialise GQL and GCE2
u0_lh   = randn(ComplexF64,2)
u0_lc   = [u0_lh[1], u0_lh[2]^2]

@show u0_lh, u0_lc

p       = [-1.0,0.5,0.5]
tspan   = (0.0,1000.0)

# solve GQL analogue
prob_lh = ODEProblem(lh!,u0_lh,tspan,p)
sol_lh  = solve(prob_lh,RK4(),adaptive=false,dt=0.005)

# solve GCE2 analogue
prob_lc = ODEProblem(lc!,u0_lc,tspan,p)
sol_lc  = solve(prob_lc,RK4(),adaptive=false,dt=0.005)

# second cumulant and power
f(x,y)  = (x,y^2)
g(x,y)  = (x,abs(y))
h(x,y)  = (x,abs(y^2))


pyplot()

t1 = sol_lh.t[:]
t2 = sol_lc.t[:]

l1 = [abs(x) for x in sol_lh[1,:]]
l2 = [abs(x) for x in sol_lc[1,:]]

c1 = [abs(y^2) for y in sol_lh[2,:]]
c2 = [abs(z) for z in sol_lc[2,:]]

dl = norm(l1 - l2)
dc = norm(c1 - c2)

@show dl,dc

# pl_lh   = Plots.plot(t1,l1,linewidth=1,xaxis="t",yguide="abs(x)",title="low: GQL",legend=false)
# pl_lc   = Plots.plot(t2,l2,linewidth=1,xaxis="t",yguide="abs(x)",title="low: GCE2",legend=false)
# pc_lh   = Plots.plot(t1,c1,linewidth=1,xaxis="t",yguide="abs(y^2)",title="cumulant: GQL",legend=false)
# pc_lc   = Plots.plot(t2,c2,linewidth=1,xaxis="t",yguide="abs(z)",title="cumulant: GCE2",legend=false)
#
# Plots.plot(Plots.plot(pl_lh,pl_lc,layout=(1,2)),Plots.plot(pc_lh,pc_lc,layout=(1,2)),layout=(2,1))

pl_lh   = Plots.plot(sol_lh,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GQL",legend=false)
pl_lc   = Plots.plot(sol_lc,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GCE2",legend=false)
pc_lh   = Plots.plot(sol_lh,vars=(h,0,2),linewidth=1,xaxis="t",yaxis="abs(y^2)",title="cumulant: GQL",legend=false)
pc_lc   = Plots.plot(sol_lc,vars=(g,0,2),linewidth=1,xaxis="t",yaxis="abs(z)",title="cumulant: GCE2",legend=false)

Plots.plot(Plots.plot(pl_lh,pl_lc,layout=(1,2)),Plots.plot(pc_lh,pc_lc,layout=(1,2)),layout=(2,1))
