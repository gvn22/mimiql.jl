using DifferentialEquations
using Plots, PyPlot

# 2D GQL analogue: low mode & high mode
function lh!(du,u,p,t)

    x,y   = u
    α,β,γ = p

    du[1] = dx = α*x + β*(x*x + y*y) + γ
    du[2] = dy = α*y + β*(x*y)

end

# 2D GCE2 analogue: low mode & second cumulant
function lc!(du,u,p,t)

    x,z   = u
    α,β,γ = p

    du[1] = dx = α*x + β*(x*x + z) + γ
    du[2] = dz = 2.0*(α*z + β*x*z)

end

# initialise GQL and GCE2
u0_gql  = randn(ComplexF64,2)
u0_gce2[1] = u0_gql[1]
u0_gce2[2] = u0_gql[2]^2

p       = [1.0,-2.0,1.5]
tspan   = (0.0,100.0)

# solve GQL equations
prob    = ODEProblem(lh!,u0_gql,tspan,p)
gql     = solve(prob,RK4())

# solve GCE2 equations
prob    = ODEProblem(lc!,u0_gce2,tspan,p)
gce2    = solve(prob,RK4())

# second cumulant, and power
f(x,y)  = (x,y*y)
g(x,y)  = (x,Float64(y*conj(y)))
h(x,y)  = (x,Float64(y*y*conj(y*y)))

pyplot()

pl_gql  = plot(gql,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="x",title="low: GQL",legend=false)
pl_gce2 = plot(gce2,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="x",title="low: GCE2",legend=false)
pc_gql  = plot(gql,vars=(h,0,2),linewidth=1,xaxis="t",yaxis="x",title="cumulant: GQL",legend=false)
pc_gce2 = plot(gce2,vars=(g,0,2),linewidth=1,xaxis="t",yaxis="y",title="cumulant: GCE2",legend=false)

plot(plot(pl_gql,pl_gce2,layout=(1,2)),plot(pc_gql,pc_gce2,layout=(1,2)),layout=(2,1))
