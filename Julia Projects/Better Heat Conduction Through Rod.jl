using Revise

using Unitful, DifferentialEquations

using Plots

#using ModelingToolkit
#using DomainSets

function rod_ODE_system!(dT_dt, T, p, t)
    density = p.density
    heat_gen = p.heat_gen
    k = p.k
    Cp = p.Cp
    walls_heat_trans_coeff = p.walls_heat_trans_coeff
    T_amb = p.T_amb
    perimeter = p.perimeter
    cross_sect_area = p.cross_sect_area
    
    N = length(T)
    dT_dt[1] = 0.0u"K/s"
    dT_dt[N] = 0.0u"K/s"

    for i in 2:(N-1)
        #println((k / (Cp * density)))
        #println((T[i+1] - 2*T[i] + 0u"K") / delta_z^2)
        #println((heat_gen - walls_heat_trans_coeff * perimeter * (T[i] - T_amb)) / (Cp * density * cross_sect_area))
        dT_dt[i] = (k / (Cp * density)) * (T[i+1] - 2*T[i] + 0u"K") / (delta_z^2) + (heat_gen - walls_heat_trans_coeff * perimeter * (T[i] - T_amb)) / (Cp * density * cross_sect_area)
    end
end

t_min = 0.0u"s"
t_max = 300.0u"s"

t_span = (t_min, t_max)

perimeter_val = 0.1u"m"
params = (
    density = 1000u"kg/m^3",
    heat_gen = 1000u"W/m^1",
    k = 400u"W/(m*K)",
    Cp = 4.0u"J/(g*K)",
    walls_heat_trans_coeff = 0.1u"W/(m^2*K)",
    T_amb = 300u"K",
    perimeter = perimeter_val,
    cross_sect_area = (perimeter_val / (2*π))^2 * π
)

z_min = 0u"m"
z_max = 1u"m"
z_steps = 50

z_vals = collect(z_min:+(z_max / z_steps):z_max)
T_vals = zeros(length(z_vals)) .* 1.0u"K" .+ 300u"K"

delta_z = (z_max - z_min) / (z_steps - 1)

rod_ODE_problem = ODEProblem(rod_ODE_system!, T_vals, t_span, params)

rod_sol = DifferentialEquations.solve(rod_ODE_problem, Tsit5())

stripped_times = ustrip.(rod_sol.t)

plt1 = plot(z_vals, rod_sol.u[1], xlabel="Position (m)", ylabel="Temperature (K)", label="T(z, t=$(round(stripped_times[1], digits=2))s)")

for i in 2:length(rod_sol.u)
    plt1 = plot!(z_vals, rod_sol.u[i], label=false)
    #label="T(z, t=$(round(stripped_times[i], digits=2))s)")
end

display(plt1)