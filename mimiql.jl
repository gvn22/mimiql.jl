using DifferentialEquations, LinearAlgebra, ParameterizedFunctions
using Plots, PlotlyJS
using Sundials

# lh -> 2D GQL analogue with one low & one high mode
# lc -> 2D GCE2 analogue with one low mode & second cumulant

lh = @ode_def_bare GQL begin

    dx = α + β*(x*x + y*y) - x
    dy = β*(2.0*x*y) - y

end α β

lc = @ode_def_bare GCE2 begin

    dx = α + β*(x*x + z) - x
    dz = 2.0*(-z + β*(2.0*x*z))

end α β

# solve over region in parameter space
αs = [10.0^i for i=-2:1:0]
βs = [10.0^i for i=-2:1:1]
lα,lβ = length(αs),length(βs)
ldiffs,cdiffs = zeros(Float64,lα,lβ),zeros(Float64,lα,lβ)

tspan = (0.0,100.0)

println(" ") # flush

for i ∈ CartesianIndices((1:lα,1:lβ))

    println("Loading parameters: α = ",αs[i[1]], " β = ",βs[i[2]])

    S = 0
    while S < 100

        u0_lh   = rand(ComplexF64,2)
        u0_lc   = [u0_lh[1], u0_lh[2]^2]

        # println("IC: x = ",u0_lh[1], " y = ",u0_lh[2])

        p       = [αs[i[1]],βs[i[2]]]

        prob_lh = ODEProblem(lh,u0_lh,tspan,p)
        sol_lh  = solve(prob_lh,RK4(),dense=true,abstol=1e-8,reltol=1e-8)

        prob_lc = ODEProblem(lc,u0_lc,tspan,p)
        sol_lc  = solve(prob_lc,RK4(),dense=true,abstol=1e-8,reltol=1e-8)

        normdx = norm([abs(sol_lh(j)[1]) - abs(sol_lc(j)[1]) for j=0:1.0:tspan[2]])
        normdz = norm([abs(sol_lh(j)[2]^2) - abs(sol_lc(j)[2]) for j=0:1.0:tspan[2]])

        if(!isnan(normdx) && !isnan(normdz))

            ldiffs[i] += normdx
            cdiffs[i] += normdz

            S += 1

        end

        # println("Diff: low = ", dl, " cumulant = ",dc)

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

    ldiffs[i] /= S
    cdiffs[i] /= S

end

pyplot()

ys = [string(i) for i = αs]
xs = [string(i) for i = βs]

p1 = wireframe(xs,ys,lnorms,zaxis="x",xaxis="β",yaxis="α",title="RK4")
p2 = wireframe(xs,ys,cnorms,zaxis="z",xaxis="β",yaxis="α")
Plots.plot(p1,p2,layout=(1,2))
