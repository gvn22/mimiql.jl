using DifferentialEquations
using Plots

function mimic2d!(du,u,p,t)

x,y = u
α,β,σ,γ,τ,δ = p

du[1] = dx = α*x + σ*(x*x + y*y) + τ
du[2] = dy = β*y + γ*x*y + 1.0 - δ

end

u0 = [0.005,0.000002]
p = [0.05,-0.25,-0.001,0.002,0.005,1.5]
tspan = (0.0,1000.0)
prob = ODEProblem(mimic2d!,u0,tspan,p)

sol = solve(prob,Tsit5(),adaptive=true)

f(x,y) = (x,y^2)

px = plot(sol,vars=1,linewidth=1,xaxis="t",yaxis="x",legend=false)
py = plot(sol,vars=2,linewidth=1,xaxis="t",yaxis="y",legend=false)
pz = plot(sol,vars=(f,0,2),linewidth=1,xaxis="t",yaxis="z",legend=false)
plot(plot(px,py,layout=(1,2)),plot(pz),layout=(2,1))
