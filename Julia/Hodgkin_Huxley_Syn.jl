using Plots, DifferentialEquations

## Function Definitions

α_n(dV) = @. 0.032*(-50 - dV)/(exp((-50 - dV)/5)- 1)
β_n(dV) = @. 0.5/(exp((-55 - dV)/40))
α_m(dV) = @. 0.32*(-52 - dV)/(exp((-52 - dV)/4) - 1)
β_m(dV) = @. 0.28*(25 + dV)/(exp((25 + dV)/5) - 1)
α_h(dV) = @. 0.128/(exp((-48 - dV)/18))
β_h(dV) = @. 4/(exp((-25 -dV)/5)+ 1)
"""
α_n(dV) = @. (0.1 - 0.01*dV)/(exp(1 - 0.1*dV)- 1)
β_n(dV) = @. 0.125/(exp(0.0125*dV))
α_m(dV) = @. (2.5 - 0.1*dV)/(exp(2.5 - 0.1*dV)- 1)
β_m(dV) = @. 4/(exp(dV/18))
α_h(dV) = @. 0.07/(exp(0.05*dV))
β_h(dV) = @. 1/(exp(3 - 0.1*dV)+ 1)

"""
n_∞(dV) = @. α_n(dV)/(α_n(dV)+β_n(dV))
m_∞(dV) = @. α_m(dV)/(α_m(dV)+β_m(dV))
h_∞(dV) = @. α_h(dV)/(α_h(dV)+β_h(dV))

# Heaviside step function (smooth version)
k = 10000
𝚯(dV) = @. 1/(1 + exp(-2*k*dV))

## Constants
const C_m     = 0.143       # nF, membrane capacitance
const G_NaMax = 7.15      # μS, maximum conductivity of Na channel
const G_KMax  = 1.43      # μS, maximum conductivity of K channel
const G_L     = 0.02672     # μS, leak conductivity

const V_Na =  50        # mV, Nernst voltage for Na
const V_K  = -95        # mV, Nernst voltage for K
const V_L  = -63.563    # mV, Nernst voltage for leak

const V_th = -20        # mV , threshold potential for transmitter
const S_max = 0.045     # maximal fraction of postsynaptically bound neurotransmitters
const τ_syn = 50        # ms, synaptic time scale
const κ = 0.5           # relative rate of transmitter binding and unbinding
const V_rev = -80     # postsynaptic reversal potential (approx 0 mV)

I_stim = 0.08    #nA, stimulus current

const V_r  = -65     # mV, resting potential
# Injected Current Function
I_inj(t) = (0 < t) & (t < 100000) ? I_stim : 0

function HH_model(du,u,p,t)
    n, m, h,S1, R1, Vm = u

    V_diff = Vm #- V_r          # difference between the rest voltage and current membrane voltage
    G_K  = G_KMax  * n^4       # Sodium conductance
    G_Na = G_NaMax * h * m^3   # Potasium conductance

    # Update transfer rate coefficients, n, m, and h
    du[1] = α_n(V_diff)*(1-n) - β_n(V_diff)*n
    du[2] = α_m(V_diff)*(1-m) - β_m(V_diff)*m
    du[3] = α_h(V_diff)*(1-h) - β_h(V_diff)*h

    #Synaptic Current | neuron 1
    du[4] = ((1/τ_syn)*(R1 - κ*S1)*(S_max - S1)/S_max)
    du[5] = (1/τ_syn)*(𝚯(Vm - V_th) - R1)

    # Update cell membrane voltage, Vm
    du[6] = (I_inj(t)  + (V_Na - Vm)*G_Na + (V_K - Vm)*G_K + (V_L - Vm)*G_L)/C_m
end

## Run Model:
u0 = [n_∞(0); m_∞(0) ; h_∞(0);0.;0.; 0]
tspan = (0.0,2000.0)
prob = ODEProblem(HH_model, u0, tspan)
sol = solve(prob,Tsit5(), saveat=0.1)

## Plotting
p1 = plot(sol.t, sol[6,:], legend=false, lw=2, ylabel="Voltage [mV]")
p2 = plot(sol.t, sol[4,:], legend=false, lw=2, ylabel="Voltage [mV]")
#p2 = plot(sol.t, I_inj.(sol.t), legend=false, lc=:red, lw=2, ylabel="Current")
p3 = plot(sol.t, sol[1:3,:]', label=["n" "m" "h"], legend=:topright, lw=2,
        xlabel="Time [ms]", ylabel="Fraction Active")

#l = grid(3, 1, heights=[0.4, 0.2 ,0.4])
#plot(p1, p2, p3, layout = l, size=(800,500))
