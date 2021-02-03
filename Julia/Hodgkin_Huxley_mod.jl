using DifferentialEquations
using Plots
using SciPy
using ProgressBars

## Function Definitions
α_n(dV) = @. (0.1 - 0.01*dV)/(exp(1 - 0.1*dV)- 1)
β_n(dV) = @. 0.125/(exp(0.0125*dV))
α_m(dV) = @. (2.5 - 0.1*dV)/(exp(2.5 - 0.1*dV)- 1)
β_m(dV) = @. 4/(exp(dV/18))
α_h(dV) = @. 0.07/(exp(0.05*dV))
β_h(dV) = @. 1/(exp(3 - 0.1*dV)+ 1)

n_∞(dV) = @. α_n(dV)/(α_n(dV)+β_n(dV))
m_∞(dV) = @. α_m(dV)/(α_m(dV)+β_m(dV))
h_∞(dV) = @. α_h(dV)/(α_h(dV)+β_h(dV))

# Heaviside step function (smooth version)
k = 10000
𝚯(dV) = @. 1/(1 + exp(-2*k*dV))

## Constants
const C_m     = 1       # μF/cm², membrane capacitance
const G_NaMax = 120     # mS/cm², maximum conductivity of Na channel
const G_KMax  = 36      # mS/cm², maximum conductivity of K channel
const G_L     = 0.3     # mS/cm², leak conductivity
const V_r  = -65        # mV, resting potential
const V_Na =  50        # mV, Nernst voltage for Na
const V_K  = -77        # mV, Nernst voltage for K
const V_L  = -54.387    # mV, Nernst voltage for leak

const V_th = -20        # mV , threshold potential for transmitter
const S_max = 0.045     # maximal fraction of postsynaptically bound neurotransmitters
const τ_syn = 0.1     # 100 μs, synaptic time scale
const κ = 0.5           # relative rate of transmitter binding and unbinding
const V_rev = -80     # postsynaptic reversal potential (approx 0 mV)

const I_stim = 10

# Injected Current Function
I_inj(t) = (0 < t) & (t < 20000) ? I_stim : 0

## Hodgkin Huxley 3 neurons model
function HH_model(du,u,p,t)
    n1, m1, h1,  S1, R1, n2, m2, h2, S2, R2, n3, m3, h3,  S3, R3, Vm1, Vm2, Vm3 = u

    V_diff1 = Vm1 - V_r # difference between the rest voltage and current membrane voltage
    V_diff2 = Vm2 - V_r
    V_diff3 = Vm3 - V_r

    G_K1  = G_KMax  * n1^4    # Sodium conductance
    G_K2  = G_KMax  * n2^4
    G_K3  = G_KMax  * n3^4

    G_Na1 = G_NaMax * h1 * m1^3   # Potasium conductance
    G_Na2 = G_NaMax * h2 * m2^3
    G_Na3 = G_NaMax * h3 * m3^3

    G12 = G21 = 30 # 30 nS --> 0.03 μS Synaptic Conductance
    G13 = G31 = 30
    G23 = G32 = 30

    # Update transfer rate coefficients, n, m, and h
    du[1] = α_n(V_diff1)*(1-n1) - β_n(V_diff1)*n1
    du[2] = α_m(V_diff1)*(1-m1) - β_m(V_diff1)*m1
    du[3] = α_h(V_diff1)*(1-h1) - β_h(V_diff1)*h1

    #Synaptic Current | neuron 1
    du[4] = ((1/τ_syn)*(R1 - κ*S1)*(S_max - S1)/S_max)
    du[5] = (1/τ_syn)*(𝚯(Vm1 - V_th) - R1)
    I_syn1 = G12*R2*(Vm1 - V_rev) + G13*R3*(Vm1 - V_rev)

    # Update transfer rate coefficients, n, m, and h | neuron 2
    du[6] = α_n(V_diff2)*(1-n2) - β_n(V_diff2)*n2
    du[7] = α_m(V_diff2)*(1-m2) - β_m(V_diff2)*m2
    du[8] = α_h(V_diff2)*(1-h2) - β_h(V_diff2)*h2

    #Synaptic Current | neuron 2
    du[9] = ((1/τ_syn)*(R2 - κ*S2)*(S_max - S2)/S_max)
    du[10] = (1/τ_syn)*(𝚯(Vm2 - V_th) - R2)
    I_syn2 =  G23*R3*(Vm2 - V_rev) + G21*R1*(Vm2 - V_rev)

    # Update transfer rate coefficients, n, m, and h | neuron 3
    du[11] = α_n(V_diff3)*(1-n3) - β_n(V_diff3)*n3
    du[12] = α_m(V_diff3)*(1-m3) - β_m(V_diff3)*m3
    du[13] = α_h(V_diff3)*(1-h3) - β_h(V_diff3)*h3

    #Synaptic Current | neuron 3
    du[14] = ((1/τ_syn)*(R3 - κ*S3)*(S_max - S3)/S_max)
    du[15] = (1/τ_syn)*(𝚯(Vm3 - V_th) - R3)
    I_syn3 = G31*R1*(Vm3 - V_rev) + G32*R2*(Vm3 - V_rev)

    # Update cell membrane voltage, Vm
    du[16] = (I_inj(t) - I_syn1 + (V_Na - Vm1)*G_Na1 + (V_K - Vm1)*G_K1 + (V_L - Vm1)*G_L)/C_m
    du[17] = (I_inj(t) - I_syn2 + (V_Na - Vm2)*G_Na2 + (V_K - Vm2)*G_K2 + (V_L - Vm2)*G_L)/C_m
    du[18] = (I_inj(t) - I_syn3 + (V_Na - Vm3)*G_Na3 + (V_K - Vm3)*G_K3 + (V_L - Vm3)*G_L)/C_m

end

## Run Model:
u0 = [n_∞(0); m_∞(0) ; h_∞(0);  0.; 0.; n_∞(0); m_∞(0) ; h_∞(0);  0.; 0.; n_∞(0); m_∞(0) ; h_∞(0);  0.; 0.; -65.1; -65.1; -65.1]
tspan = (0.0,200.0)
prob = ODEProblem(HH_model, u0, tspan)
sol = solve(prob, saveat = 0.01, reltol = 1e-10,abstol = 1e-8, progress = true)

## Plotting
p1 = plot(sol.t, sol[5,:], legend=false, lw=2, ylabel="Voltage [mV]")
p2 = plot!(sol.t, sol[10,:], legend=false, lw=2, ylabel="Voltage [mV]")
p3 = plot!(sol.t, sol[15,:], legend=false, lw=2, ylabel="Voltage [mV]")
#pall = plot(sol[5,:], sol[10,:], sol[15,:])
#savefig("3Rate_WLC.png")

p_1 = plot(sol.t, sol[16,:], legend=false, lw=2, ylabel="Voltage [mV]")
p_2 = plot!(sol.t, sol[17,:], legend=false, lw=2, ylabel="Voltage [mV]")
p_3 = plot!(sol.t, sol[18,:], legend=false, lw=2, ylabel="Voltage [mV]")
p_all = plot(sol[16,:], sol[17,:], sol[18,:])
#plot(p_1,p_2,p_3, layout = (3,1))
#savefig("3HH_WLC.png")

#p2 = plot(sol.t, I_inj.(sol.t), legend=false, lc=:red, lw=2, ylabel="Current")
#p3 = plot(sol.t, sol[1:3,:]', label=["n" "m" "h"], legend=:topright, lw=2, xlabel="Time [ms]", ylabel="Fraction Active")

## Frequency-Current curves

### Find Peaks
function find_peak(series)
    # residence time vector for each neuron (LV equation)
    peak_u, _ = SciPy.signal.find_peaks(series, height = 20, distance = 5)
    #width_peak_u, _ = SciPy.signal.peak_widths(series,peak_u, rel_height = 0.8)
    return peak_u
end

function find_peak_R(series)
    # residence time vector for each neuron (LV equation)
    peak_u, _ = SciPy.signal.find_peaks(series, height = 0.8, distance = 5)
    #width_peak_u, _ = SciPy.signal.peak_widths(series,peak_u, rel_height = 0.8)
    return peak_u
end

function freq(peak_series)
    freq_Hz = 1/((peak_series[4] - peak_series[3])*0.01*0.001)
    return freq_Hz
end


#peaks_HH1_10 = find_peak(sol[16,:])
#freq_HH1_10 = freq(peaks_HH1_10) #Frequency in Hz
#peaks_HH2_10 = find_peak(sol[17,:])
#peaks_HH3_10 = find_peak(sol[18,:])

## looping over different I_inj values
I_values = collect(10:50)
freqs_HH1 = []
freqs_R1 = []

for I_v in tqdm(I_values, total = length(I_values))
    # Injected Current Function
    I_inj(t) = (0 < t) & (t < 2000) ? I_v : 0

    function HH_model(du,u,p,t)
        n1, m1, h1,  S1, R1, n2, m2, h2, S2, R2, n3, m3, h3,  S3, R3, Vm1, Vm2, Vm3 = u

        V_diff1 = Vm1 - V_r # difference between the rest voltage and current membrane voltage
        V_diff2 = Vm2 - V_r
        V_diff3 = Vm3 - V_r

        G_K1  = G_KMax  * n1^4    # Sodium conductance
        G_K2  = G_KMax  * n2^4
        G_K3  = G_KMax  * n3^4

        G_Na1 = G_NaMax * h1 * m1^3   # Potasium conductance
        G_Na2 = G_NaMax * h2 * m2^3
        G_Na3 = G_NaMax * h3 * m3^3

        G12 = G21 = 30 # 30 nS --> 0.03 μS Synaptic Conductance
        G13 = G31 = 30
        G23 = G32 = 30

        # Update transfer rate coefficients, n, m, and h
        du[1] = α_n(V_diff1)*(1-n1) - β_n(V_diff1)*n1
        du[2] = α_m(V_diff1)*(1-m1) - β_m(V_diff1)*m1
        du[3] = α_h(V_diff1)*(1-h1) - β_h(V_diff1)*h1

        #Synaptic Current | neuron 1
        du[4] = ((1/τ_syn)*(R1 - κ*S1)*(S_max - S1)/S_max)
        du[5] = (1/τ_syn)*(𝚯(Vm1 - V_th) - R1)
        I_syn1 = G12*R2*(Vm1 - V_rev) + G13*R3*(Vm1 - V_rev)

        # Update transfer rate coefficients, n, m, and h | neuron 2
        du[6] = α_n(V_diff2)*(1-n2) - β_n(V_diff2)*n2
        du[7] = α_m(V_diff2)*(1-m2) - β_m(V_diff2)*m2
        du[8] = α_h(V_diff2)*(1-h2) - β_h(V_diff2)*h2

        #Synaptic Current | neuron 2
        du[9] = ((1/τ_syn)*(R2 - κ*S2)*(S_max - S2)/S_max)
        du[10] = (1/τ_syn)*(𝚯(Vm2 - V_th) - R2)
        I_syn2 =  G23*R3*(Vm2 - V_rev) + G21*R1*(Vm2 - V_rev)

        # Update transfer rate coefficients, n, m, and h | neuron 3
        du[11] = α_n(V_diff3)*(1-n3) - β_n(V_diff3)*n3
        du[12] = α_m(V_diff3)*(1-m3) - β_m(V_diff3)*m3
        du[13] = α_h(V_diff3)*(1-h3) - β_h(V_diff3)*h3

        #Synaptic Current | neuron 3
        du[14] = ((1/τ_syn)*(R3 - κ*S3)*(S_max - S3)/S_max)
        du[15] = (1/τ_syn)*(𝚯(Vm3 - V_th) - R3)
        I_syn3 = G31*R1*(Vm3 - V_rev) + G32*R2*(Vm3 - V_rev)

        # Update cell membrane voltage, Vm
        du[16] = (I_inj(t) - I_syn1 + (V_Na - Vm1)*G_Na1 + (V_K - Vm1)*G_K1 + (V_L - Vm1)*G_L)/C_m
        du[17] = (I_inj(t) - I_syn2 + (V_Na - Vm2)*G_Na2 + (V_K - Vm2)*G_K2 + (V_L - Vm2)*G_L)/C_m
        du[18] = (I_inj(t) - I_syn3 + (V_Na - Vm3)*G_Na3 + (V_K - Vm3)*G_K3 + (V_L - Vm3)*G_L)/C_m
        #println(I_syn1)
    end

    u0 = [n_∞(0); m_∞(0) ; h_∞(0);  0.; 0.; n_∞(0); m_∞(0) ; h_∞(0);  0.; 0.; n_∞(0); m_∞(0) ; h_∞(0);  0.; 0.; -65.1; -65.1; -65.1]
    tspan = (0.0,200.0)
    prob = ODEProblem(HH_model, u0, tspan)
    sol = solve(prob, saveat = 0.01, reltol = 1e-10,abstol = 1e-8, progress = true)

    peak1 = find_peak(sol[16,:])
    freq1 = freq(peak1)

    peakR1 = find_peak_R(sol[5,:])
    freqR1 = freq(peakR1)

    append!(freqs_HH1, freq1)
    append!(freqs_R1, freqR1)
end

"""
plot(I_values, freqs_HH1)
savefig("I_vs_freq_HH1.pdf")
plot(I_values, freqs_R1)
savefig("I_vs_freq_R1.pdf")
"""

## Export files for better plots

"""
open("Neuron1HH_90.txt", "w") do f
  for i in sol[16,:]
    println(f, i)
  end
end #

open("Neuron2HH_90.txt", "w") do f
  for i in sol[17,:]
    println(f, i)
  end
end #

open("Neuron3HH_90.txt", "w") do f
  for i in sol[18,:]
    println(f, i)
  end
end #
"""
