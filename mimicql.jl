using DifferentialEquations
using Plots, PyPlot

# gql analogue with only low and high modes
function lh!(du,u,p,t)

x,y = u
α,γ = p

du[1] = dx = α*x + (x*x + y*y) + γ
du[2] = dy = α*y + x*y

end

u0 = randn(ComplexF64,2)
# p = [0.05,-0.25,-0.001,0.002,0.0005,-0.000005]
p = [-0.05,0.005]
tspan = (0.0,100.0)

prob = ODEProblem(lh!,u0,tspan,p)
sol = solve(prob,RK4(),adaptive=true)

# second cumulant analogue in 2d
f(x,y) = (x,Float64(y*conj(y)))

pyplot()
px = plot(sol,vars=1,linewidth=1,xaxis="t",yaxis="x",legend=false)
py = plot(sol,vars=2,linewidth=1,xaxis="t",yaxis="y",legend=false)
pz = plot(sol,vars=(f,0,2),linewidth=1,xaxis="t",yaxis="z = y^2",legend=false)
plot(plot(px,py,layout=(1,2)),plot(pz),layout=(2,1))

# gql analogue with a low and two high modes
# function lhh!(du,u,p,t)
#
# x,y,z = u
# α,β,σ,γ,τ,δ = p
#
# β2,γ2,δ2 = β/2.0,γ/3.0,δ/4.0
#
# du[1] = dx = α*x + σ*(x*x + y*y + y*z + z*z) + τ
# du[2] = dy = β*y + γ*(x*y + x*z) + δ
# du[3] = dz = β2*z + γ2*(x*z + x*y) + δ2
#
# end

# u0 = [0.005,0.000005,0.001]
# p = [0.05,-0.25,-0.001,0.002,0.0005,-0.000005]
# tspan = (0.0,1000.0)
# prob = ODEProblem(lhh!,u0,tspan,p)
# sol = solve(prob,Tsit5(),adaptive=true)

# second cumulant analogue in 3d
# f(x,y,z) = (x,y^2 + y*z + z^2)
# px = plot(sol,vars=1,linewidth=1.5,xaxis="t",yaxis="x",legend=false)
# py = plot(sol,vars=2,linewidth=1.5,xaxis="t",yaxis="y",legend=false)
# pz = plot(sol,vars=3,linewidth=1.5,xaxis="t",yaxis="z",legend=false)
# pyz = plot(sol,vars=(f,0,2,3),linewidth=1.5,xaxis="t",yaxis="z = y^2",legend=false)
# plot(px,py,pz,layout=(3,1))
