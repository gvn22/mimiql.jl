using DifferentialEquations
using Plots

function mimic!(du,u,p,t)

x,y = u
α,β,γ,μ,τ = p

x0 = 0.005

du[1] = α*x - y*y + x0
du[2] = α*x*y - β*y

end

u0 = [1.0005,0.2]
p = [0.05,0.0025,0.5,2.0,0.5]
tspan = (0.0,100.0)
prob = ODEProblem(mimic!,u0,tspan,p)

sol = solve(prob,Tsit5())

# plot(sol,linewidth=1,layout=(2,1),legend=false)

gr()
f(x,y) = (x,y^2)

px = plot(sol,vars=1,linewidth=1,xaxis="t",yaxis="x",legend=false)
py = plot(sol,vars=2,linewidth=1,xaxis="t",yaxis="y",legend=false)
pz = plot(sol,vars=(f,0,2),linewidth=1,xaxis="t",yaxis="z",legend=false)
plot(plot(px,py,layout=(1,2)),plot(pz),layout=(2,1))
