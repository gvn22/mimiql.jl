using DifferentialEquations, LinearAlgebra, ParameterizedFunctions, OrdinaryDiffEq
using Plots, PlotlyJS

# f(xy) -> 2D GQL analogue with one low & one high mode
# g(xz) -> 2D GCE2 analogue with one low mode & second cumulant

fxy = @ode_def_bare begin

    dx = a - x + c*(x*x + y*y)
    dy = - y + c*(2.0*x*y)

end a c

gxz = @ode_def_bare begin

    dx = a - x + c*(x*x + z)
    dz = 2.0*(-z + c*(2.0*x*z))

end a c

function plot_at_param(pa,pb)

    u0_xy   = rand(ComplexF64,2)
    u0_xz   = [u0_xy[1], u0_xy[2]^2]
    p       = [pa,pb]
    tspan   = (0.0,100.0)

    prob_xy = ODEProblem(fxy,u0_xy,tspan,p)
    sol_xy  = solve(prob_xy,RK4(),dense=true,abstol=1e-8,reltol=1e-8)

    prob_xz = ODEProblem(gxz,u0_xz,tspan,p)
    sol_xz  = solve(prob_xz,RK4(),dense=true,abstol=1e-8,reltol=1e-8)

    normdx = norm([abs(sol_xy(j)[1]) - abs(sol_xz(j)[1]) for j=0:1.0:tspan[2]])
    normdz = norm([abs(sol_xy(j)[2]^2) - abs(sol_xz(j)[2]) for j=0:1.0:tspan[2]])

    @show normdx,normdz

    g(x,y)  = (x,abs(y))
    h(x,y)  = (x,abs(y^2))

    plotly()

    pl_xy   = Plots.plot(sol_xy,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GQL",legend=false)
    pl_xz   = Plots.plot(sol_xz,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GCE2",legend=false)
    pc_xy   = Plots.plot(sol_xy,vars=(h,0,2),linewidth=1,xaxis="t",yaxis="abs(y^2)",title="cumulant: GQL",legend=false)
    pc_xz   = Plots.plot(sol_xz,vars=(g,0,2),linewidth=1,xaxis="t",yaxis="abs(z)",title="cumulant: GCE2",legend=false)

    pcomp   = Plots.plot(Plots.plot(pl_xy,pl_xz,layout=(1,2)),Plots.plot(pc_xy,pc_xz,layout=(1,2)),layout=(2,1))

    Plots.display(pcomp)

end

# solve over region in parameter space
as = [10.0^i for i=-2:1:1]
cs = [10.0^i for i=-2:1:1]
la,lc = length(as),length(cs)
xdiffs,zdiffs = zeros(Float64,la,lc),zeros(Float64,la,lc)

println(" ") # flush output

for i ∈ CartesianIndices((1:la,1:lc))

    a,c = as[i[1]],cs[i[2]]
    println("Loading parameters: a = ", a, " c = ", c)

    S = 0.0

    while S < 100.0

        u0_xy   = randn(ComplexF64,2)
        u0_xz   = [u0_xy[1], u0_xy[2]^2]
        tspan   = (0.0,100.0)
        p       = [a,c]

        prob_xy = ODEProblem(fxy,u0_xy,tspan,p)
        sol_xy  = solve(prob_xy,RK4(),dense=true,abstol=1e-8,reltol=1e-8)

        prob_xz = ODEProblem(gxz,u0_xz,tspan,p)
        sol_xz  = solve(prob_xz,RK4(),dense=true,abstol=1e-8,reltol=1e-8)

        normdx = norm([abs(sol_xy(j)[1]) - abs(sol_xz(j)[1]) for j=0:1.0:tspan[2]])
        normdz = norm([abs(sol_xy(j)[2]^2) - abs(sol_xz(j)[2]) for j=0:1.0:tspan[2]])

        if(!isnan(normdx) && !isnan(normdz))

            xdiffs[i] += normdx
            zdiffs[i] += normdz

            S += 1.0

        end

    end

    xdiffs[i] /= S
    zdiffs[i] /= S

    println("Difference: dx = ", xdiffs[i], " dz = ", zdiffs[i])

end

plotly()

xs = [string(i) for i ∈ as]
ys = [string(i) for i ∈ cs]

p1 = Plots.plot(xs,ys,xdiffs',st=:contourf,color=:matter,xaxis="a",yaxis="c",title="Δx")
# p2 = Plots.plot(xs,ys,zdiffs',st=:contourf,color=:matter,xaxis="a",yaxis="c",title="Δz")
# Plots.plot(p1,p2,layout=(1,2))

# plot_at_param(0.1, 10.0)


# using Distributed
# using DifferentialEquations
# using Plots
#
# addprocs(2)
# @everywhere using DifferentialEquations
# @everywhere using ParameterizedFunctions,LinearAlgebra
#
# @everywhere function prob_func(prob,i,repeat)
#   remake(prob,u0=rand()*prob.u0)
# end
#
# g(x,y)  = (x,abs(y))
#
# u0 = rand(ComplexF64,2)
# # Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0
# prob = ODEProblem(gxz,u0,(0.0,100.0),[1.0,10.0])
# ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
# sim = solve(ensemble_prob,RK4(),EnsembleThreads(),trajectories=100)
# pyplot()
# Plots.plot(sim,vars=(g,0,1),linealpha=0.4)
# summ = EnsembleSummary(sim,vars=(g,0,1),0:0.1:10)
# pyplot() # Note that plotly does not support ribbon plots
# Plots.plot(summ,fillalpha=0.5)
