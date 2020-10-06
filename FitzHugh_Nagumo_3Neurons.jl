using DifferentialEquations
using Plots; pyplot()
using MultivariateStats


function fitzhugh_nagumo(du,u,p,t)
  ### 6D fitzhugh_nagumo system
  x_A, y_A,  x_B, y_B,  x_C, y_C= u
  α, β, γ, δ = p

  du[1] = dx_A = α*(x_A - (x_A*x_A*x_A)/3 + y_A)+ β*(x_A - x_B)
  du[2] = dy_A = -(1/α)*(x_A + δ*y_A - γ)


  du[3] = dx_B = α*(x_B - (x_B*x_B*x_B)/3  + y_B) + β*(x_B - x_C)
  du[4] = dy_B = -(1/α)*(x_B + δ*y_B - γ)

  du[5] = dx_C = α*(x_C - (x_C*x_C*x_C)/3 + y_C) + β*(x_C  - x_A)
  du[6] = dy_C = -(1/α)*(x_C + δ*y_C - γ)

end

# Defining initial conditions and time
u0 = [1.5,1.2, 1.0, 1.3, 1.1, 1.6]
tspan = (0.0,10.0)
dt = 0.1
p = [0.7,1.5,0.7,1.5]

## ODE problem, take a look at DifferentialEquation.jl docs
prob_1 = ODEProblem(fitzhugh_nagumo,u0,tspan,p, saveat =dt )

#solution
sol = solve(prob_1)
u_s = Array(sol)


#plotting solutions
plot(sol)
##### Performing PCA
M = fit(PCA, u_s; maxoutdim = 3)

result = transform(M, u_s)


#savefig("solutions.png")

# PCA plot
plot(result[1,:], result[2,:], result[3,:])
#plot(u_s[1,:], u_s[3,:], u_s[5,:])
#savefig("phase_portrait.png")

dres1 = vcat([result[1, i+1] - result[1,i] for i in 1:size(u_s)[2] -1], 0)
dres2 = vcat([result[2, i+1] - result[2,i] for i in 1:size(u_s)[2] -1], 0)
dres3 = vcat([result[3, i+1] - result[3,i] for i in 1:size(u_s)[2] -1], 0)

"""
anim = @animate for i ∈ 1:size(u_s)[2]
    Plots.quiver(result[1,i:end], result[2,i:end], result[3,i:end], quiver = (dres1[i:end], dres2[i:end], dres3[i:end]))
end
gif(anim, "WLC_FHN.gif", fps = 15)
"""
