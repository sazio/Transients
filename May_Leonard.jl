using DifferentialEquations
using Plots; pyplot()


function may_leonard(du,u,p,t)
  ### 3D lotka volterra system
  ## http://www.math.nthu.edu.tw/~sbhsu/May-Leonard.pdf
  x, y, z= u
  α_1, β_1, α_2, β_2, α_3, β_3 = p
  du[1] = dx = x*(1 - x - α_1*y - β_1*z)
  du[2] = dy = y*(1 - y - β_2*x - α_2*z)
  du[3] = dz = z*(1 - z - α_3*x - β_3*y)
end

# Defining initial conditions and time
u0 = [1.5,1.2, 1.0]
tspan = (0.0,1000.0)
dt = 0.1
p = [2.0,0.5,2.0,0.5,2.0,0.5]

## ODE problem, take a look at DifferentialEquation.jl docs
prob = ODEProblem(may_leonard,u0,tspan,p, saveat =dt )

#solution
sol = solve(prob)
u_s = Array(sol)

du_s1 = vcat([u_s[1, i+1] - u_s[1,i] for i in 1:size(u_s)[2] -1], 0)
du_s2 = vcat([u_s[2, i+1] - u_s[2,i] for i in 1:size(u_s)[2] -1], 0)
du_s3 = vcat([u_s[3, i+1] - u_s[3,i] for i in 1:size(u_s)[2] -1], 0)
#my_cgrad = cgrad([:red, :yellow, :blue])
#plotting
plot(sol)
savefig("solutions.png")
plot(u_s[1,:], u_s[2,:], u_s[3,:])
savefig("phase_portrait.png")



#surf(u_s[1,:], u_s[2,:], u_s[3,:], supp = [du_s1 du_s2 du_s3], w = "vectors filled head", lw = 2)
anim = @animate for i ∈ 1:size(u_s)[2]
    Plots.quiver(u_s[1,i:end], u_s[2,i:end], u_s[3,i:end], quiver = (du_s1[i:end], du_s2[i:end], du_s3[i:end]))
end
gif(anim, "WLC_long.gif", fps = 15)

# Phase Portrait
#Plots.quiver(u_s[1,:], u_s[2,:], u_s[3,:], quiver = (du_s1, du_s2, du_s3))
