using DifferentialEquations
using Plots;
using MultivariateStats
using SciPy
using ProgressBars

function Total_LV(du,u,p,t, r = 2.2)
  ### Total Compeition Lotka Volterra System --> WLC
  ### In this case we have a map  !!!!!
  ## http://www.math.nthu.edu.tw/~sbhsu/May-Leonard.pdf
  q, v, w, x, y, z = u

  ρ_12, ρ_13, ρ_14, ρ_15, ρ_16, ρ_21, ρ_23, ρ_24, ρ_25, ρ_26, ρ_31, ρ_32, ρ_34, ρ_35, ρ_36, ρ_41, ρ_42, ρ_43, ρ_45, ρ_46, ρ_51, ρ_52, ρ_53, ρ_54, ρ_56, ρ_61, ρ_62, ρ_63, ρ_64, ρ_65 = p

   du[1] = dq = r*q*(1 - q - ρ_12*v - ρ_13*w - ρ_14*x - ρ_15*y - ρ_16*z)
   du[2] = dv = r*v*(1 - ρ_21*q -v -ρ_23*w - ρ_24*x - ρ_25*y - ρ_26*z)
   du[3] = dw  = r*w*(1 - ρ_31*q - ρ_32*v -w - ρ_34*x - ρ_35*y - ρ_36*z)
   du[4] = dx = r*x*(1 - ρ_41*q - ρ_42*v - ρ_43*w - x - ρ_45*y - ρ_46*z)
   du[5] = dy = r*y*(1 - ρ_51*q - ρ_52*v - ρ_53*w - ρ_54*x - y - ρ_56*z)
   du[6] = dz = r*z*(1 - ρ_61*q - ρ_62*v - ρ_63*w - ρ_64*x - ρ_65*y - z)

end




# Defining initial conditions and time
u0 = [0.0001, 0.0001, 0.0001, 0.0001, 0.04, 0.0001]
tspan = (0.0,100000.0)
dt = 1
p = [1.05, 1.001, 1.001, 1.001, 0, 0, 1.15, 1.001, 1.001, 1.001, 1.001, 0, 1.14, 1.001, 1.001, 1.001, 1.001, 0, 1.07, 1.001, 1.001, 1.001, 1.001, 0, 1.09, 1.001, 1.001, 1.001, 1.001, 0]

## ODE problem, take a look at DifferentialEquation.jl docs
prob = DiscreteProblem(Total_LV,u0,tspan,p, saveat =dt)
#solution
sol = solve(prob)
u_s = Array(sol)
#plotting
plot(sol, layout = (3,2))
#savefig("solutions.png")


##### Performing PCA for phase portrait
#=
pca_fit = fit(PCA, u_s; maxoutdim = 3)

result = transform(pca_fit, u_s)


# PCA plot
plot(result[1,:], result[2,:], result[3,:])

savefig("total_LV_6eqs.png")
=#

###### Residence time vectors


function res_time_vec(series)
    # residence time vector for each neuron (LV equation)
    peak_u, _ = SciPy.signal.find_peaks(series, height = 0.0008, distance = 60)
    width_peak_u, _ = SciPy.signal.peak_widths(series,peak_u, rel_height = 0.8)

    return width_peak_u
end


res_time_1 = Int.(floor.(res_time_vec(u_s[1,:])))
res_time_2 = Int.(floor.(res_time_vec(u_s[2,:])))
res_time_3 = Int.(floor.(res_time_vec(u_s[3,:])))
res_time_4 = Int.(floor.(res_time_vec(u_s[4,:])))
res_time_5 = Int.(floor.(res_time_vec(u_s[5,:])))
res_time_6 = Int.(floor.(res_time_vec(u_s[6,:])))
res_time = [res_time_1, res_time_2, res_time_3, res_time_4, res_time_5, res_time_6]

plot(res_time, layout = (3,2))

D_ = 6

function diff_res_time(res_time, k = 1 ,  D = D_)

    diff = res_time[k + D] - res_time[k]
    return diff
end


# Scheme 1
function Scheme1(res_time, D = D_ )

    scheme = Vector{Int64}()
    for i in 1:Int(floor(size(res_time_1)[1]/D) -1)
        append!(scheme, diff_res_time(res_time, D*i))
    end
    return scheme
end


# Scheme 2
function Scheme2(res_time, D = D_ )

    scheme = Vector{Int64}()
    for i in 1:(size(res_time)[1] - D)
        append!(scheme, diff_res_time(res_time, i))
    end
    return scheme
end

#plot(Scheme2(res_time_3))



power_spectrum_3 = SciPy.signal.periodogram(Scheme2(res_time_3))
plot(power_spectrum_3)
#---------------
r_s = 1:0.001:3

spectra_3 = power_spectrum_3
for r_i in ProgressBar(r_s)

    function Total_LV(du,u,p,t, r = r_i)
      ### Total Compeition Lotka Volterra System --> WLC
      ### In this case we have a map  !!!!!
      ## http://www.math.nthu.edu.tw/~sbhsu/May-Leonard.pdf
      q, v, w, x, y, z = u

      ρ_12, ρ_13, ρ_14, ρ_15, ρ_16, ρ_21, ρ_23, ρ_24, ρ_25, ρ_26, ρ_31, ρ_32, ρ_34, ρ_35, ρ_36, ρ_41, ρ_42, ρ_43, ρ_45, ρ_46, ρ_51, ρ_52, ρ_53, ρ_54, ρ_56, ρ_61, ρ_62, ρ_63, ρ_64, ρ_65 = p

       du[1] = dq = r*q*(1 - q - ρ_12*v - ρ_13*w - ρ_14*x - ρ_15*y - ρ_16*z)
       du[2] = dv = r*v*(1 - ρ_21*q -v -ρ_23*w - ρ_24*x - ρ_25*y - ρ_26*z)
       du[3] = dw  = r*w*(1 - ρ_31*q - ρ_32*v -w - ρ_34*x - ρ_35*y - ρ_36*z)
       du[4] = dx = r*x*(1 - ρ_41*q - ρ_42*v - ρ_43*w - x - ρ_45*y - ρ_46*z)
       du[5] = dy = r*y*(1 - ρ_51*q - ρ_52*v - ρ_53*w - ρ_54*x - y - ρ_56*z)
       du[6] = dz = r*z*(1 - ρ_61*q - ρ_62*v - ρ_63*w - ρ_64*x - ρ_65*y - z)

    end

    prob = DiscreteProblem(Total_LV,u0,tspan,p, saveat =dt)
    #solution
    sol = solve(prob)
    u_s = Array(sol)
    # Residence time
    res_time_3 = Int.(floor.(res_time_vec(u_s[3,:])))
    # Power Spectrum
    power_spectrum_3 = SciPy.signal.periodogram(Scheme2(res_time_3))
    if r_i == 1.0
        spectra_3 = power_spectrum_3
    else
        spectra_3 = vcat(spectra_3, power_spectrum_3)
    end

end


a_s = r_s[32:end]
spectra_final = [s for s in spectra_3 if !all(isempty.(s))]
x_s = [ a for a  in a_s]

plot(x_s,first.(spectra_final), last.(spectra_final), st=:surface)

size(a_s)
