using DifferentialEquations, LinearAlgebra, ParameterizedFunctions
using Plots, PlotlyJS

# lh -> 2D GQL analogue with one low & one high mode
# lc -> 2D GCE2 analogue with one low mode & second cumulant

lh = @ode_def_bare begin

    dx = α + β*(x*x + y*y) - x
    dy = β*(2.0*x*y) - y

end β α

lc = @ode_def_bare begin

    dx = α + β*(x*x + z) - x
    dz = 2.0*(-z + β*(2.0*x*z))

end β α

# solve over a region in parameter space
# βs,αs   = [0.1,1.0,10.0],[0.01,0.1,1.0]
βs = [10.0^i for i=-2:1:1]
αs = [10.0^i for i=-2:1:1]

px,py   = length(βs),length(αs)

lnorms  = zeros(Float64,px,py)
cnorms  = zeros(Float64,px,py)

S = 100

for i ∈ CartesianIndices((1:px,1:py))

    println(" ")
    println("Params: β = ",βs[i[1]], " α = ",αs[i[2]])

    tspan   = (0.0,100.0)

    for j=1:1:S

        u0_lh   = randn(ComplexF64,2)
        u0_lc   = [u0_lh[1], u0_lh[2]^2]

        # println("IC: x = ",u0_lh[1], " y = ",u0_lh[2])

        p       = [βs[i[1]],αs[i[2]]]

        prob_lh = ODEProblem(lh,u0_lh,tspan,p)
        sol_lh  = solve(prob_lh,Tsit5(),adaptive=false,dt=0.001,abstol=1e-6,reltol=1e-6)

        prob_lc = ODEProblem(lc,u0_lc,tspan,p)
        sol_lc  = solve(prob_lc,Tsit5(),adaptive=false,dt=0.001,abstol=1e-6,reltol=1e-6)

        l1 = [abs(x) for x in sol_lh[1,:]]
        l2 = [abs(x) for x in sol_lc[1,:]]

        c1 = [abs(y^2) for y in sol_lh[2,:]]
        c2 = [abs(z) for z in sol_lc[2,:]]

        dl = norm(l1 - l2)
        dc = norm(c1 - c2)
        if(!isnan(dl) && !isnan(dc))
            lnorms[i] += dl
            cnorms[i] += dc
        end

        println("Diff: low = ", dl, " cumulant = ",dc)

        # g(x,y)  = (x,abs(y))
        # h(x,y)  = (x,abs(y^2))

        # plotly()
        #
        # pl_lh   = Plots.plot(sol_lh,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GQL",legend=false)
        # pl_lc   = Plots.plot(sol_lc,vars=(g,0,1),linewidth=1,xaxis="t",yaxis="abs(x)",title="low: GCE2",legend=false)
        # pc_lh   = Plots.plot(sol_lh,vars=(h,0,2),linewidth=1,xaxis="t",yaxis="abs(y^2)",title="cumulant: GQL",legend=false)
        # pc_lc   = Plots.plot(sol_lc,vars=(g,0,2),linewidth=1,xaxis="t",yaxis="abs(z)",title="cumulant: GCE2",legend=false)
        #
        # pcomp   = Plots.plot(Plots.plot(pl_lh,pl_lc,layout=(1,2)),Plots.plot(pc_lh,pc_lc,layout=(1,2)),layout=(2,1))
        #
        # Plots.display(pcomp)

    end

    lnorms[i] /= S
    cnorms[i] /= S

end

pyplot()

xs = [string(i) for i = βs]
ys = [string(i) for i = αs]

p1 = wireframe(xs,ys,lnorms',zaxis="x",xaxis="β",yaxis="α",title="Tsit5")
p2 = wireframe(xs,ys,cnorms',zaxis="z",xaxis="β",yaxis="α")
# Plots.plot(p1,p2,layout=(1,2))
