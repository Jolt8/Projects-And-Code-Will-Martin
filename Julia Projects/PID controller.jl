using Revise
using Plots
#using DifferentialEquations
#using NLsolve

SP = 400 #K
u_bias = 0.0

total_time = 100

time_step = 5

K_c = 8 #proportional gain

T_i = 10 #Integral time

T_D = 0.1 #Derivative time

u_bias = 0

K_p = 2

function gen_PV(PV_old, u_input, time_step)
    process_gain = 2
    process_time_constant = 100
    return PV_old + time_step * ((1 / process_time_constant) * (-PV_old + process_gain * u_input))
end

num_steps = (total_time / time_step)

num_steps = Int(total_time / time_step)
PV = zeros(num_steps + 1)
e = zeros(num_steps + 1)
u = zeros(num_steps + 1)

PV[1] = 300

integral_sum = 0
last_error = 0

for t in 1:num_steps
    e[t] = SP - PV[t]

    P_out = K_c * e[t]

    integral_sum += e[t] * time_step
    I_out = (K_c / T_i) * integral_sum

    D_out = K_c * T_D * ((e[t] - last_error) / time_step)

    u[t] = P_out + I_out + D_out + u_bias

    PV[t+1] = gen_PV(PV[t], u[t], time_step)
    last_error = e[t]
end

u[end] = u[end-1]

for t in eachindex(u)
    println("Time (t): ", t * time_step, ", Error (e): ", e[t], ", Process Variable, ", PV[t], ", Controller Output (u): ", u[t])
end

ts = []
for t in eachindex(u)
    push!(ts, t)
end

plt = plot(ts, e,
    xlabel="Time",
    ylabel="Process Variable",
    legend=:best, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    label="e",
    framestyle=:box)

plot!(ts, PV, 
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    label="PV",)

plot!(ts, u,
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    label="u",)



