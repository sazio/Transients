#Pkg.add("DifferentialEquations")
using DifferentialEquations
using Plots #pyplot()
using MultivariateStats


function step_G(x)
  if x <= 0
    return 0.0
  else
    return 1.0
  end
end

function fitzhugh_nagumo(du,u,p,t)
  ### 9D fitzhugh_nagumo system
  x_A, y_A,z_A,  x_B, y_B, z_B, x_C, y_C, z_C= u
  τ_1, ν_min, S, a, b, τ_2 = p

  du[1] = dx_A = (1/τ_1)*(x_A - (x_A*x_A*x_A)/3 - y_A + -z_A*(x_A - ν_min) + 0.35 +S)
  du[2] = dy_A = (x_A - b*y_A + a)
  du[3] = dz_A = (1/τ_2)*( 2*step_G(x_C) - z_A)

  du[4] = dx_B = (1/τ_1)*(x_B - (x_B*x_B*x_B)/3 - y_B + -z_B*(x_B - ν_min) + 0.35 +S)
  du[5] = dy_B = (x_B - b*y_B + a)
  du[6] = dz_B = (1/τ_2)*(2*step_G(x_A)  -z_B)

  du[7] = dx_C = (1/τ_1)*(x_C - (x_C*x_C*x_C)/3 - y_C + -z_C*(x_C - ν_min) + 0.35 +S)
  du[8] = dy_C = (x_C - b*y_C + a)
  du[9] = dz_C = (1/τ_2)*(2*step_G(x_B) -z_C)
end

# Defining initial conditions and time
#u0 = [1.1, 0.78, 0.43, -1.34, 0.55, -0.27, 1.23, -0.67, 0.21]
#resting state
u0 = [-1.4,-0.52, 0.2,-1.35,-0.73, 0.3,-1.3, -0.72, 0.1]

tspan = (0.0, 200.0)
dt = 0.01
p = [0.08,-1.5,0.15,0.7, 0.8, 3.1]

## ODE problem, take a look at DifferentialEquation.jl docs
prob_1 = ODEProblem(fitzhugh_nagumo,u0,tspan,p, saveat =dt)

#solution
sol = solve(prob_1, reltol = 1e-10)
u_s = Array(sol)


#plotting solutions
#plot(sol)
#savefig("solutions.png")

##### Performing PCA

pca_fit = fit(PCA, u_s; maxoutdim = 3)
result = transform(pca_fit, u_s)

# PCA plot
plot(result[1,:], result[2,:], result[3,:])
#plot(result1, result2, result3)
#plot(u_s[1,:], u_s[3,:], u_s[5,:])
savefig("phase_portraitFHN.pdf")

"""
dres1 = vcat([result[1, i+1] - result[1,i] for i in 1:size(u_s)[2] -1], 0)
dres2 = vcat([result[2, i+1] - result[2,i] for i in 1:size(u_s)[2] -1], 0)
dres3 = vcat([result[3, i+1] - result[3,i] for i in 1:size(u_s)[2] -1], 0)


anim = @animate for i ∈ 1:size(u_s)[2]
    Plots.quiver(result[1,i:end], result[2,i:end], result[3,i:end], quiver = (dres1[i:end], dres2[i:end], dres3[i:end]))
end
gif(anim, "WLC_FHN.gif", fps = 15)

# Phase Portrait
Plots.quiver(result[1,:], result[2,:], result[3,:], quiver = (dres1, dres2, dres3))
"""
