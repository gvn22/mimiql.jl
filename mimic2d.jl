using DifferentialEquations
using Plots

function mimic!(du,u,p,t)

x,y,z = u
σ,ρ,β = p
du[1] = dx = σ*(y-x)
du[2] = dy = x*(ρ-z) - y
du[3] = dz = x*y - β*z

end

u0 = [1.0,2.0,1.0]
p = [10.0,28.0,8.0/3.0]
tspan = (0.0,100.0)
prob = ODEProblem(mimic!,u0,tspan,p)

sol = solve(prob)

plot(sol,linewidth=1,xaxis="t",layout=(3,1),legend=false)
