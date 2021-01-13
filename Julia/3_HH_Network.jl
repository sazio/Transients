using Plots, DifferentialEquations

## Function Definitions
α_n(dV) = @. 0.032*(50 - dV)/(exp((-5 - dV)/5)- 1)
β_n(dV) = @. 0.5/(exp((-55 - dV)/40))

α_m(dV) = @. 0.32*(-52 - dV)/(exp((-52 - dV)/4) - 1)
β_m(dV) = @. 0.28*(25 + dV)/(exp((25 + dV)/5) - 1)

α_h(dV) = @. 0.128/(exp((-48 - dV)/18))
β_h(dV) = @. 4/(exp((-25 -dV)/5)+ 1)

n_∞(dV) = @. α_n(dV)/(α_n(dV)+β_n(dV))
m_∞(dV) = @. α_m(dV)/(α_m(dV)+β_m(dV))
h_∞(dV) = @. α_h(dV)/(α_h(dV)+β_h(dV))

# Heaviside step function (smooth version)
const k = 10000
𝚯(dV) = @. 1/(1 + exp(-2*k*dV))

## Constants
const C_m     = 0.03       # μF/cm², membrane capacitance
const G_NaMax = 360     # mS/cm², maximum conductivity of Na channel
const G_KMax  = 70      # mS/cm², maximum conductivity of K channel
const G_L     = 1.     # mS/cm², leak conductivity
#const V_r  = -65        # mV, resting potential
const V_Na =  50        # mV, Nernst voltage for Na
const V_K  = -95        # mV, Nernst voltage for K
const V_L  = -64    # mV, Nernst voltage for leak

const V_th = -20        #mV , threshold potential for transmitter
const S_max = 0.045     # maximal fraction of postsynaptically bound neurotransmitters
const τ_syn = 0.05      # 50 ms
const κ = 0.5           # relative rate of transmitter binding and unbinding
const V_rev = 61     # postsynaptic reversal potential (approx 0 mV)

# Injected Current Function
I_inj(t) = (5 < t) & (t < 30) ? 10 : 0

function HH_model(du,u,p,t)
    n1, m1, h1,  S1, R1, n2, m2, h2, S2, R2, n3, m3, h3,  S3, R3, Vm1, Vm2, Vm3 = u

    V_diff1 = Vm1 #- V_r # difference between the rest voltage and current membrane voltage
    V_diff2 = Vm2 #- V_r
    V_diff3 = Vm3 #- V_r

    G_K1  = G_KMax  * n1^4    # Sodium conductance
    G_K2  = G_KMax  * n2^4
    G_K3  = G_KMax  * n3^4

    G_Na1 = G_NaMax * h1 * m1^3   # Potasium conductance
    G_Na2 = G_NaMax * h2 * m2^3
    G_Na3 = G_NaMax * h3 * m3^3

    G12 = G21 = 1 # Synpatic Conductance
    G13 = G31 = 1
    G23 = G32 = 1

    # Update transfer rate coefficients, n, m, and h | neuron 1
    du[1] = α_n(V_diff1)*(1-n1) - β_n(V_diff1)*n1
    du[2] = α_m(V_diff1)*(1-m1) - β_m(V_diff1)*m1
    du[3] = α_h(V_diff1)*(1-h1) - β_h(V_diff1)*h1

    #Synaptic Current | neuron 1
    I_syn1 = G12*S1*(Vm2 - V_rev) + G13*S1*(Vm3 - V_rev)
    du[4] = (1/τ_syn)*(R1 - κ*S1)*(S_max - S1)/S_max
    du[5] = (1/τ_syn)*(𝚯(Vm1 - V_th) - R1)

    # Update transfer rate coefficients, n, m, and h | neuron 2
    du[6] = α_n(V_diff2)*(1-n2) - β_n(V_diff2)*n2
    du[7] = α_m(V_diff2)*(1-m2) - β_m(V_diff2)*m2
    du[8] = α_h(V_diff2)*(1-h2) - β_h(V_diff2)*h2

    #Synaptic Current | neuron 2
    I_syn2 = G21*S2*(Vm1 - V_rev) + G23*S2*(Vm3 - V_rev)
    du[9] = (1/τ_syn)*(R2 - κ*S2)*(S_max - S2)/S_max
    du[10] = (1/τ_syn)*(𝚯(Vm2 - V_th) - R2)

    # Update transfer rate coefficients, n, m, and h | neuron 3
    du[11] = α_n(V_diff3)*(1-n3) - β_n(V_diff3)*n3
    du[12] = α_m(V_diff3)*(1-m3) - β_m(V_diff3)*m3
    du[13] = α_h(V_diff3)*(1-h3) - β_h(V_diff3)*h3

    #Synaptic Current | neuron 3
    I_syn3 = G31*S3*(Vm1 - V_rev) + G32*S3*(Vm2 - V_rev)
    du[14] = (1/τ_syn)*(R3 - κ*S3)*(S_max - S3)/S_max
    du[15] = (1/τ_syn)*(𝚯(Vm3 - V_th) - R3)

    # Update cell membrane voltage, Vm
    du[16] = (I_inj(t) + I_syn1  + (V_Na - Vm1)*G_Na1 + (V_K - Vm1)*G_K1 + (V_L - Vm1)*G_L)/C_m
    du[17] = (I_inj(t) + I_syn2 + (V_Na - Vm2)*G_Na2 + (V_K - Vm2)*G_K2 + (V_L - Vm2)*G_L)/C_m
    du[18] = (I_inj(t) + I_syn3 + (V_Na - Vm3)*G_Na3 + (V_K - Vm3)*G_K3 + (V_L - Vm3)*G_L)/C_m
end

## Run Model:
u0 = [n_∞(0); m_∞(0) ; h_∞(0);  0.5; 0.5; n_∞(0); m_∞(0) ; h_∞(0);  0.5; 0.5; n_∞(0); m_∞(0) ; h_∞(0);  0.5; 0.5; -65.9; -65.9; -65.9]
tspan = (0.0,50.0)
prob = ODEProblem(HH_model, u0, tspan)
sol = solve(prob, saveat=0.1)

sol[16,:]


## Plotting
p1 = plot(sol.t, sol[16,:], legend=false, lw=2, ylabel="Voltage [mV]")
p_2 = plot(sol.t, sol[17,:], legend=false, lw=2, ylabel="Voltage [mV]")
p_3 = plot(sol.t, sol[18,:], legend=false, lw=2, ylabel="Voltage [mV]")




p2 = plot(sol.t, I_inj.(sol.t), legend=false, lc=:red, lw=2, ylabel="Current")
p3 = plot(sol.t, sol[1:3,:]', label=["n" "m" "h"], legend=:topright, lw=2,
        xlabel="Time [ms]", ylabel="Fraction Active")

l = grid(3, 1, heights=[0.4, 0.2 ,0.4])
plot(p1, p2, p3, layout = l, size=(800,500))
