using NeuralPDE
using Quadrature, Cubature, Cuba
using Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Plots
using PyPlot
#GalacticOptim
# 3D PDE
@parameters x y t θ
@variables u(..)
@derivatives Dx'~x
@derivatives Dy'~y
@derivatives Dt'~t

α = 2.0
β = 0.5


# 3D PDE - Solution of Liouville's theorem for May-Leonard
eq  = Dt(u(x,y,t,θ)) +  Dx(u(x,y,t,θ))*((α + β)/2)*x*(1-x-2*y) + Dy(u(x,y,t,θ))*((α + β)/2)*y*(1-2*x-y) ~ 0


# Initial and boundary conditions

A0   = 1.0
muX0 = 0.
muY0 = 0.
sgX0 = 1.0
sgY0 = 1.0

bcs = [u(x,y,0,θ) ~ A0 * exp(-(((x-muX0)^2/(2*sgX0^2))) +((y-muY0)^2/(2*sgY0^2)))]

# Space and time domains
domains = [x ∈ IntervalDomain(0.0,1.0),
           y ∈ IntervalDomain(0.0,1.0),
           t ∈ IntervalDomain(0.0,10.0)]

# Discretization
dx = 0.1; dy = 0.1; dt = 0.1

# Neural network
chain = FastChain(FastDense(3,16,Flux.σ),FastDense(16,16,Flux.σ),FastDense(16,1))

# Quadrature training -- new, by Kirill
q_strategy = NeuralPDE.QuadratureTraining(algorithm =CubaCuhre(),reltol=1e-8,abstol=1e-8,maxiters=100)
discretization = NeuralPDE.PhysicsInformedNN([dx,dy,dt],chain,strategy = q_strategy)


pde_system = PDESystem(eq,bcs,domains,[x,y,t],[u])
prob = discretize(pde_system, discretization)

cb = function (p,l)
    println("Current loss is: $l")
    return false
end

res = GalacticOptim.solve(prob, GalacticOptim.ADAM(0.03), progress = true, cb = cb, maxiters=4000)
phi = discretization.phi



####

xs,ys,ts = [domain.domain.lower:dx:domain.domain.upper for domain in domains]

u_predict = [reshape([first(phi([x,y,t],res.minimizer)) for x in xs for y in ys], (length(xs),length(ys))) for t in ts]


#analytic_sol_func(x,t) = [exp(-t)*sin(pi*x), exp(-t)*cos(pi*x), (1+pi^2)*exp(-t)]
#u_real  = [reshape([first(analytic_sol_func(x,t)[i] for x in xs for y in ys], (length(xs),length(ys))) for t in ts]
#u_real  = [[analytic_sol_func(x,t)[i] for t in ts for x in xs] for i in 1:3]


#u_predict = [reshape([first(phi([x,t],res.minimizer)) for x in xs for t in ts])]


maxlim = maximum(maximum(u_predict[t]) for t = 1:length(ts))
minlim = minimum(minimum(u_predict[t]) for t = 1:length(ts))


result = @animate for time = 1:length(ts)
    Plots.plot(xs, ys, u_predict[time],st=:surface,camera=(30,30), zlim=(minlim,maxlim), clim=(minlim,maxlim),
                title = string("ψ: max = ",round(maxlim, digits = 3)," min = ", round(minlim, digits = 3),"\\n"," t = ",time))
end

gif(result, "Liouville.gif", fps = 6)

result = @animate for time = 1:length(ts)
    Plots.plot(xs, ys, u_predict[time],st=:contour,camera=(30,30), zlim=(minlim,maxlim), clim=(minlim,maxlim),
                title = string("ψ: max = ",round(maxlim, digits = 3)," min = ", round(minlim, digits = 3),"\\n"," t = ",time))
end

gif(result, "liouville_contour.gif", fps = 6)
