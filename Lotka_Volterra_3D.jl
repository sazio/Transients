using DifferentialEquations
using Plots; pyplot()


function lotka_volterra(du,u,p,t)
  ### 3D lotka volterra system
  x, y, z= u
  α, β, δ, γ, f, g, e = p
  du[1] = dx = α*x - β*x*y
  du[2] = dy = -δ*y + γ*x*y - e*y*z
  du[3] = dz = -f*z + g*y*z
end

# Defining initial conditions and time
u0 = [1.0,1.0, 1.0]
tspan = (0.0,10.0)
dt = 0.0001
p = [1.5,1.0,2.0,1.0,1.50,2.0,2.0]

## ODE problem, take a look at DifferentialEquation.jl docs
prob = ODEProblem(lotka_volterra,u0,tspan,p, saveat =dt )

#solution
sol = solve(prob)
u_s = Array(sol)

#my_cgrad = cgrad([:red, :yellow, :blue])
#plotting
plot(sol)
plot(u_s[1,:], u_s[2,:], u_s[3,:])
