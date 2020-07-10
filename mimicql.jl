using DifferentialEquations
using Plots, PyPlot

# 2D GQL analogue: low mode & high mode
function lh!(du,u,p,t)

    x,y   = u
    α,β,γ = p

    du[1] = dx = α*x + β*(x*x + y*y) + γ*(x - 1.0)
    du[2] = dy = α*y + β*(2.0*x*y)

end

# 2D GCE2 analogue: low mode & second cumulant
function lc!(du,u,p,t)

    x,z   = u
    α,β,γ = p

    du[1] = dx = α*x + β*(x*x + z) + γ*(x - 1.0)
    du[2] = dz = 2.0*(α*z + β*(2.0*x*z))

end

# initialise GQL and GCE2
u0_lh   = randn(ComplexF64,2)
u0_lc   = [u0_lh[1], u0_lh[2]^2]

# @show u0_lh, u0_lc

p       = [-1.0,1.5,0.25]
tspan   = (0.0,100.0)

# solve GQL analogue
prob_lh = ODEProblem(lh!,u0_lh,tspan,p)
sol_lh  = solve(prob_lh,Tsit5(),reltol=1e-8)

# solve GCE2 analogue
prob_lc = ODEProblem(lc!,u0_lc,tspan,p)
sol_lc  = solve(prob_lc,Tsit5(),reltol=1e-8)

# second cumulant and power
f(x,y)  = (x,y^2)
g(x,y)  = (x,abs(y)^2)
h(x,y)  = (x,abs(y^2)^2)

pyplot()

pl_lh   = plot(sol_lh,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="x",title="low: GQL",legend=false)
pl_lc   = plot(sol_lc,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="x",title="low: GCE2",legend=false)
pc_lh   = plot(sol_lh,vars=(h,0,2),linewidth=1,xaxis="t",yaxis="y^2",title="cumulant: GQL",legend=false)
pc_lc   = plot(sol_lc,vars=(g,0,2),linewidth=1,xaxis="t",yaxis="z",title="cumulant: GCE2",legend=false)

plot(plot(pl_lh,pl_lc,layout=(1,2)),plot(pc_lh,pc_lc,layout=(1,2)),layout=(2,1))
