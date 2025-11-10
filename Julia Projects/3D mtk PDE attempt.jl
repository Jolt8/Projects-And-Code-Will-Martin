using ModelingToolkit
using DifferentialEquations
using Plots
using GLMakie

import ModelingToolkit: t_nounits as t, D_nounits as D

#this was taken from the MTK standard library 
@connector HeatPort begin
    @parameters begin
        T_guess = 273.15 + 20
        Q_flow_guess = 0.0
    end

    @variables begin
        T(t), [guess = T_guess]
        Q_flow(t), [guess = Q_flow_guess, connect = Flow]
    end
end

# A simple 0D heat capacitor, representing one node.
@mtkmodel HeatCapacitor begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        C = 1.0, [description = "Heat capacity of element"]
    end
    @variables begin
        T(t) = 293.15, [description = "Temperature of element [K]"]
        der_T(t) = 0.0, [description = "Time derivative of temperature [K/s]"]
    end

    @equations begin
        T ~ port.T
        D(T) ~ port.Q_flow / C
    end
end

@mtkmodel Element1D begin
    @components begin
        port_a = HeatPort()
        port_b = HeatPort()
    end
    @variables begin
        dT(t), [guess = 0.0]
        Q_flow(t), [guess = 0.0]
    end
    @equations begin
        dT ~ port_a.T - port_b.T
        port_a.Q_flow ~ Q_flow
        port_a.Q_flow + port_b.Q_flow ~ 0
    end
end

# A 0D thermal conductor, representing the connection between nodes.
@mtkmodel ThermalConductor begin
    @extend Q_flow, dT = element1d = Element1D()
    @parameters begin
        G = 1.0, [description = "Thermal conductance [W/K]"]
    end
    @equations begin
        Q_flow ~ G * dT
    end
end

# A constant temperature boundary condition
@mtkmodel SomeFixedTemperature begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        T_fixed = 303.15, [description = "Fixed temperature [K]"] # 30 °C
    end
    @equations begin
        port.T ~ T_fixed
    end
end

# A constant heat flow boundary condition
@mtkmodel SomeHeatSource begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        Q_flow_fixed = 10.0, [description = "Fixed heat flow [W]"]
    end
    @equations begin
        port.Q_flow ~ Q_flow_fixed
    end
end

rows = 3
cols = 3
layers = 3
n_nodes = rows * cols * layers

nodes = [HeatCapacitor(name=Symbol("node_", i, "_", j, "_", k), C=1/n_nodes) for i in 1:rows, j in 1:cols, k in 1:layers]

v_conds = [ThermalConductor(name=Symbol("v_cond_", i, "_", j, "_", k)) for i in 1:(rows-1), j in 1:cols, k in 1:layers]
h_conds = [ThermalConductor(name=Symbol("h_cond_", i, "_", j, "_", k)) for i in 1:rows, j in 1:(cols-1), k in 1:layers]
d_conds = [ThermalConductor(name=Symbol("d_cond_", i, "_", j, "_", k)) for i in 1:rows, j in 1:cols, k in 1:(layers-1)]
#d for depth I guess

# Add boundary condition components
left_bcs = [SomeHeatSource(name=Symbol("heat_source_", i, "_", k)) for i in 1:rows, k in 1:layers]
right_bcs = [SomeFixedTemperature(T_fixed=293.15, name=Symbol("fixed_temp_", i, "_", k)) for i in 1:rows, k in 1:layers]

connections = Equation[]

[push!(connections, connect(nodes[i, j, k].port, v_conds[i, j, k].port_a)) for i in 1:(rows-1), j in 1:cols, k in 1:layers]
[push!(connections, connect(v_conds[i, j, k].port_b, nodes[i+1, j, k].port)) for i in 1:(rows-1), j in 1:cols, k in 1:layers]

[push!(connections, connect(nodes[i, j, k].port, h_conds[i, j, k].port_a)) for i in 1:rows, j in 1:(cols-1), k in 1:layers]
[push!(connections, connect(h_conds[i, j, k].port_b, nodes[i, j+1, k].port)) for i in 1:rows, j in 1:(cols-1), k in 1:layers]

[push!(connections, connect(nodes[i, j, k].port, d_conds[i, j, k].port_a)) for i in 1:rows, j in 1:cols, k in 1:(layers-1)]
[push!(connections, connect(d_conds[i, j, k].port_b, nodes[i, j, k+1].port)) for i in 1:rows, j in 1:cols, k in 1:(layers-1)]

[push!(connections, connect(left_bcs[i, k].port, nodes[i, 1, k].port)) for i in 1:rows, k in 1:layers]
[push!(connections, connect(right_bcs[i, k].port, nodes[i, cols, k].port)) for i in 1:rows, k in 1:layers]

# --- Simulation Setup ---
eqs_total = System[]
all_systems = vcat(
    vec(nodes),
    vec(h_conds),
    vec(v_conds),
    vec(d_conds),
    vec(left_bcs),
    vec(right_bcs)
)
all_systems

@named rod = ODESystem(connections, t, systems=all_systems)

sys = structural_simplify(rod)  

tspan = (0.0, 10.0)

prob = ODEProblem(sys, [], tspan)

sol = solve(prob)
"""
initial_temps = sol.u[1]
final_temps = sol.u[end]
T_min = minimum(initial_temps)
T_max = maximum(final_temps)

fig = Figure(size = (800, 700))
ax = Axis3(fig[1, 1],
    title = "Temperature Distribution in 2D Heat Plate",
    xlabel = "Column (X)",
    ylabel = "Row (Y)",
    zlabel = "Temperature (K)"
)

solution_matrix = reshape(vcat(sol.u[24], (ones(rows) * 293.15)), rows, cols)
temp_obs = Observable(solution_matrix)
heatmap = GLMakie.heatmap(1:4, 1:4, temp_obs)

for i in eachindex(sol.u)
    solution_matrix = reshape(vcat(sol.u[i], (ones(rows) * 293.15)), rows, cols)
    temp_obs = Observable(solution_matrix)
    heatmap = GLMakie.heatmap(1:4, 1:4, temp_obs)
    display(heatmap)
    sleep(0.1)
end


unknowns(rod)
vars = unknowns(rod)[1:4:(4*n_nodes)]
plot(sol,
     vars=vars,
     title="Temperature Distribution in Heat Rod",
     xlabel="Time (s)",
     ylabel="Temperature (K)",
     legend=:topleft,
     #labels=["Node " .* string.(1:10) for i in 1:10]'
    )
"""