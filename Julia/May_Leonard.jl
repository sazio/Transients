using DifferentialEquations
using Plots;


function may_leonard(du,u,p,t)
  ### 3D lotka volterra system
  ## http://www.math.nthu.edu.tw/~sbhsu/May-Leonard.pdf
  x, y, z= u
  ρ_12, ρ_13, ρ_21,ρ_23, ρ_31, ρ_32 = p
  du[1] = dx = x*(1 - x - ρ_12*y - ρ_13*z)
  du[2] = dy = y*(1 - y - ρ_21*x - ρ_23*z)
  du[3] = dz = z*(1 - z - ρ_31*x - ρ_32*y)
end

# Defining initial conditions and time
u0 = [1.0,1.5, 1.2]
tspan = (0.0,300.0)
dt = 0.1
p = [2.0,0.5,0.5,2.0,2.0,0.5]

## ODE problem, take a look at DifferentialEquation.jl docs
prob = ODEProblem(may_leonard,u0,tspan,p, saveat =dt )

#solution
sol = solve(prob)
u_s = Array(sol)

#my_cgrad = cgrad([:red, :yellow, :blue])
#plotting
plot(sol)
#savefig("solutions.png")
plot(u_s[1,:], u_s[2,:], u_s[3,:])
savefig("phase_portraitML.png")


"""
du_s1 = vcat([u_s[1, i+1] - u_s[1,i] for i in 1:size(u_s)[2] -1], 0)
du_s2 = vcat([u_s[2, i+1] - u_s[2,i] for i in 1:size(u_s)[2] -1], 0)
du_s3 = vcat([u_s[3, i+1] - u_s[3,i] for i in 1:size(u_s)[2] -1], 0)

#surf(u_s[1,:], u_s[2,:], u_s[3,:], supp = [du_s1 du_s2 du_s3], w = "vectors filled head", lw = 2)
anim = @animate for i ∈ 1:size(u_s)[2]
    Plots.quiver(u_s[1,i:end], u_s[2,i:end], u_s[3,i:end], quiver = (du_s1[i:end], du_s2[i:end], du_s3[i:end]))
end
gif(anim, "WLC_long.gif", fps = 15)

# Phase Portrait
#Plots.quiver(u_s[1,:], u_s[2,:], u_s[3,:], quiver = (du_s1, du_s2, du_s3))

"""
