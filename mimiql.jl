using DifferentialEquations, LinearAlgebra, ParameterizedFunctions, OrdinaryDiffEq
using Plots, PlotlyJS
plotly()

# f(xy) -> 2D GQL analogue with one low & one high mode
# g(xz) -> 2D GCE2 analogue with one low mode & second cumulant

fxy = @ode_def_bare begin

    dx = a - x + c*(x*x + y*y)
    dy = - y - c*(2.0*x*y)

end a c

gxz = @ode_def_bare begin

    dx = a - x + c*(x*x + z)
    dz = 2.0*(-z - c*(2.0*x*z))

end a c

# define parameter space
as = [10.0^i for i=-2:0.5:1]
cs = [10.0^i for i=-2:0.5:1]
la,lc = length(as),length(cs)
normdx,normdz = zeros(Float64,la,lc),zeros(Float64,la,lc)

println(" ") # flush output

T  = 1000.0
Δt = 0.001
Nt = Int(T/Δt) + 2

xdiffs = zeros(Float64,length(as),length(cs),Nt)
zdiffs = zeros(Float64,length(as),length(cs),Nt)

for i ∈ CartesianIndices((1:la,1:lc))

    a,c = as[i[1]],cs[i[2]]
    println("Loading parameters: a = ", a, " c = ", c)

    S  = 0.0

    while S < 1.0

        u0_xy   = randn(ComplexF64,2)
        u0_xz   = [u0_xy[1], u0_xy[2]^2]
        tspan   = (0.0,T)
        p       = [a,c]

        prob_xy = ODEProblem(fxy,u0_xy,tspan,p)
        sol_xy  = solve(prob_xy,RK4(),dt=Δt,adaptive=false,dense=false,calck=false)

        prob_xz = ODEProblem(gxz,u0_xz,tspan,p)
        sol_xz  = solve(prob_xz,RK4(),dt=Δt,adaptive=false,dense=false,calck=false)

        # normdx = norm([abs(sol_xy(j)[1]) - abs(sol_xz(j)[1]) for j=0:1.0:tspan[2]])
        # normdz = norm([abs(sol_xy(j)[2]^2) - abs(sol_xz(j)[2]) for j=0:1.0:tspan[2]])

        if(!isnan(sol_xy[1,Nt]) && !isnan(sol_xz[1,Nt]) && !isnan(sol_xy[2,Nt]) && !isnan(sol_xz[2,Nt]))

            xdiffs[i,:] = [abs(u[1]) - abs(v[1]) for (u,v) in zip(sol_xy.u,sol_xz.u)]
            zdiffs[i,:] = [abs(u[2]^2) - abs(v[2]) for (u,v) in zip(sol_xy.u,sol_xz.u)]

            g(x,y)  = (x,abs(y))
            h(x,y)  = (x,abs(y^2))

            pl_xy   = Plots.plot(sol_xy,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GQL",legend=false)
            pl_xz   = Plots.plot(sol_xz,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GCE2",legend=false)
            pc_xy   = Plots.plot(sol_xy,vars=(h,0,2),linewidth=1,xaxis="t",yaxis="abs(y^2)",title="cumulant: GQL",legend=false)
            pc_xz   = Plots.plot(sol_xz,vars=(g,0,2),linewidth=1,xaxis="t",yaxis="abs(z)",title="cumulant: GCE2",legend=false)

            pcomp   = Plots.plot(Plots.plot(pl_xy,pl_xz,layout=(1,2)),Plots.plot(pc_xy,pc_xz,layout=(1,2)),layout=(2,1),fmt=:pdf)

            fn = string("./run/sol_a",i[1],"_c",i[2])
            Plots.pdf(pcomp, fn)

            S += 1.0

        else
            println("Ran into nans. Re-running...")
        end

    end

    # xdiffs[i,:] = xdiffs[i,:]/S
    # zdiffs[i,:] = zdiffs[i,:]/S

    normdx[i] = norm(xdiffs[i,:])
    normdz[i] = norm(zdiffs[i,:])

    # println("Difference: dx = ", xdiffs[i], " dz = ", zdiffs[i])

end


xs = ["0.01","0.0316","0.1","0.316","1.0","3.16","10.0"]
ys = ["0.01","0.0316","0.1","0.316","1.0","3.16","10.0"]

p1 = Plots.plot(xs,ys,xdiffs[:,:,10000]',st=:contourf,color=:matter,xaxis="a",yaxis="c",title="Δx")
# p2 = Plots.plot(xs,ys,zdiffs[:,:,1000]',st=:contourf,color=:matter,xaxis="a",yaxis="c",title="Δz")

# function plot_at_param(pa,pb)
#
#     u0_xy   = rand(ComplexF64,2)
#     u0_xz   = [u0_xy[1], u0_xy[2]^2]
#     p       = [pa,pb]
#     tspan   = (0.0,1000.0)
#
#     prob_xy = ODEProblem(fxy,u0_xy,tspan,p)
#     sol_xy  = solve(prob_xy,RK4(),dt=0.001,adaptive=false)
#
#     prob_xz = ODEProblem(gxz,u0_xz,tspan,p)
#     sol_xz  = solve(prob_xz,RK4(),dt=0.001,adaptive=false)
#
#     normdx = norm([abs(sol_xy(j)[1]) - abs(sol_xz(j)[1]) for j=0:1.0:tspan[2]])
#     normdz = norm([abs(sol_xy(j)[2]^2) - abs(sol_xz(j)[2]) for j=0:1.0:tspan[2]])
#
#     @show normdx,normdz
#
#     g(x,y)  = (x,abs(y))
#     h(x,y)  = (x,abs(y^2))
#
#     pl_xy   = Plots.plot(sol_xy,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GQL",legend=false)
#     pl_xz   = Plots.plot(sol_xz,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GCE2",legend=false)
#     pc_xy   = Plots.plot(sol_xy,vars=(h,0,2),linewidth=1,xaxis="t",yaxis="abs(y^2)",title="cumulant: GQL",legend=false)
#     pc_xz   = Plots.plot(sol_xz,vars=(g,0,2),linewidth=1,xaxis="t",yaxis="abs(z)",title="cumulant: GCE2",legend=false)
#
#     pcomp   = Plots.plot(Plots.plot(pl_xy,pl_xz,layout=(1,2)),Plots.plot(pc_xy,pc_xz,layout=(1,2)),layout=(2,1))
#
#     Plots.display(pcomp)
#
# end

# plot_at_param(10.0, 1.0)
