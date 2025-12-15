using Unitful 
using DifferentialEquations
using Clapeyron
using Plots
using RecursiveArrayTools
using GLMakie
using SparseConnectivityTracer, ADTypes

#include("Chemeqs/src/Chemeqs.jl")
#using .Chemeqs 
#import .Chemeqs:overall_rate, arrenhius_equation_pre_exponential_factor, rate_law_for_chemicals #For some reason this only works if I do this


function arrenhius_equation_pre_exponential_factor(k, Ea, T)
    R_GAS = 8.314u"J/(mol*K)"
    result = (k / exp(-Ea / (R_GAS * T)))
    #result = ustrip(result)
    #return result * Unitful.unit(k)
    return result
end

function arrenhius_equation_rate_constant(A, Ea, T)
    R_GAS = 8.314u"J/(mol*K)"
    return (A * exp(-Ea / (R_GAS * T)))
end

struct HomogenousChemicalReaction
    name::String
    # Use Dicts for reactants/products with species name => stoichiometric coefficient
    # For rate law, these are the powers in elementary reactions.
    # For mass balance, they are the coefficients.
    reactants::Dict{String, Int}
    products::Dict{String, Int}
    #To maintain order throughout the process
    all_chemicals::Vector{String}
    kf_A::Any # Pre-exponential factor for forward reaction
    kf_Ea::Any  # Activation energy for forward reaction
    kr_A::Any # Pre-exponential factor for reverse reaction
    kr_Ea::Any   # Activation energy for reverse reaction
end
struct HeterogeneousChemicalReaction
    name::String
    # Use Dicts for reactants/products with species name => stoichiometric coefficient
    # For rate law, these are the powers in elementary reactions.
    # For mass balance, they are the coefficients.
    reactants::Dict{String, Int}
    products::Dict{String, Int}
    #To maintain order throughout the process
    all_chemicals::Vector{String}
    kf_A # Pre-exponential factor for forward reaction
    kf_Ea # Activation energy for forward reaction
    kr_A # Pre-exponential factor for reverse reaction
    kr_Ea  # Activation energy for reverse reaction
end

function K_eq(forward_k, reverse_k)
    return forward_k / reverse_k
end

#K1 is initial equilibrum constant, returns equilibrium constant after temperature change
function van_t_hoft(K_eq1, T1, T2, heat_of_reaction)
    R_GAS = 8.314u"J/(mol*K)"
    return K_eq1 * exp((-heat_of_reaction / R_GAS) * (1/T2 - 1/T1))
end
function K_gibbs_free(T_ref, T_actual, ΔG_rxn_ref, ΔH_rxn_ref)
    T_actual = uconvert(u"K", T_actual) # Ensure T is a plain number in Kelvin

    R = 8.314e-3u"kJ/(mol*K)" # Gas constant in kJ
    
    K_ref = exp(-ΔG_rxn_ref / (R * T_ref))

    ln_K_ratio = (-ΔH_rxn_ref / R) * (1/T_actual - 1/T_ref)
    
    K_T = K_ref * exp(ln_K_ratio)
    
    return ustrip(K_T) # Return a unitless equilibrium constant
end

function reversible_pressure_rate_law(reaction, temperature, concentrations_dict)
    # 1. Calculate forward and reverse rate constants at the current temperature
    kf = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
    kr = arrenhius_equation_rate_constant(reaction.kr_A, reaction.kr_Ea, temperature)

    forward_term = 1.0u"L^0" # Initialize with a unitless 1
    for (reactant, nu) in reaction.reactants
        concentration = max(concentrations_dict[reactant], 0.0u"mol/L")
        forward_term *= concentration^nu
    end

    reverse_term = 1.0u"L^0" # Initialize with a unitless 1
    for (product, nu) in reaction.products
        concentration = max(concentrations_dict[product], 0.0u"mol/L")
        reverse_term *= concentration^nu
    end
    
    net_reaction_rate = ((kf * forward_term) - (kr * reverse_term))
    #println("name: ", reaction.name, ", kf: ", kf, ", Forward Term: ", forward_term, ", kr: ", kr, ", Reverse Term: ", reverse_term, ", Net Reaction Rate: ", (kf * forward_term) - (kr * reverse_term))
    #println("")
    
    rates = Dict{String, Any}()
    for chemical in reaction.all_chemicals
        net_stoichiometry = 0
        if haskey(reaction.products, chemical)
            net_stoichiometry += reaction.products[chemical]
        end
        if haskey(reaction.reactants, chemical)
            net_stoichiometry -= reaction.reactants[chemical]
        end
        rates[chemical] = net_stoichiometry * net_reaction_rate
        #println("chemical: ", chemical, ", net_reaction_rate: ", net_reaction_rate, ", net stoichiometry: ", net_stoichiometry, ", rate for chemical: ", rates[chemical])
    end
    
    return [rates, net_reaction_rate]
end

#overall_rate(SMR_reaction, SMR_model, reactor_inlet_temperature, 100000u"Pa", SMR_initial_C)


mutable struct PIDState
    last_error::typeof(1.0u"K")
    integral_error::typeof(1.0u"K*s")
    last_time::typeof(1.0u"s")
end

mutable struct PBRModel
    reactions::Vector{HomogenousChemicalReaction}
    clapeyron_model::PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}
    components::Vector{String}
    heat_of_reactions_dict::Dict{String, typeof(1.0u"kJ/mol")}

    N_points::Int64

    pid_state::PIDState

    desired_reactor_temp::typeof(1.0u"K")

    proportional_gain::Float64
    integral_time::Float64
    derivative_time::Float64

    previous_errors::Vector{typeof(1.0u"K")}
    
    reactor_length::typeof(1.0u"m")
    reactor_diameter::typeof(1.0u"m")

    catalyst_density::typeof(1.0u"kg/m^3")
    
    ambient_overall_heat_transfer_coeff::typeof(1.0u"W/(m^2*K)")
    ambient_temperature::typeof(1.0u"K")
    
    reactor_wall_overall_heat_transfer_coeff::typeof(1.0u"W/(m^2*K)")

    jacket_fluid_cp::typeof(1.0u"J/(kg*K)")

    jacket_fluid_inlet_temperature::typeof(1.0u"K")

    heat_capacities_dict::Dict{String, typeof(1.0u"J/(mol*K)")}

    initial_F_in::Vector{typeof(1.0u"mol/s")}

    #the reaction_mixture conductivity can't be calculated with clapeyron.jl but it can with coolprop, but that isn't written in julia, 
    #sooo... we'll just use an approximation
end


function calculate_rates_and_properties(T::typeof(1.0u"K"), P, concentrations_vec, p::PBRModel)
    reactions = p.reactions
    components = p.components

    total_mol_per_l = sum(concentrations_vec)
    concentrations_dict = Dict(zip(components, concentrations_vec))
    z_vec = [concentration / total_mol_per_l for concentration in concentrations_vec]
    
    net_reaction_rates = Dict()
    total_species_rates = Dict(c => 0.0u"mol/(L*s)" for c in components)
    for i in eachindex(reactions)
        species_rates, net_rate = reversible_pressure_rate_law(reactions[i], T, concentrations_dict)

        net_reaction_rates[reactions[i].name] = net_rate

        #println(reactions[i].name, reactions[i].all_chemicals, species_rates)
        for (chemical, rate) in species_rates
            total_species_rates[chemical] += rate
        end
    end
    return (net_reaction_rates=net_reaction_rates, total_species_rates=total_species_rates, concentrations_dict=concentrations_dict, z_vec=z_vec)
end

function pbr_ode_system!(du, u, p, t)
    section_length = p.reactor_length / p.N_points

    section_volume = p.reactor_length * (pi * (p.reactor_diameter / 2)^2)

    section_wall_area = p.reactor_length * (pi * p.reactor_diameter)

    kg_catalyst_per_section = p.reactor_length * (pi * (p.reactor_diameter / 2)^2) * p.catalyst_density

    proportional_gain = p.proportional_gain
    integral_time = p.integral_time
    derivative_time = p.derivative_time

    P_vec = u.x[1]
    T_vec = u.x[2]
    T_jacket_vec = u.x[3]
    m_flow_jacket_vec = u.x[4]
    moles_mat = u.x[5]

    dP_dt = du.x[1]
    dT_dt = du.x[2]
    dT_jacket_dt = du.x[3]
    dm_flow_jacket_dt = du.x[4]
    dmoles_dt = du.x[5]

    jacket_flow_flux = fill(1.0u"W", N_points)

    for i in 1:p.N_points
        P = P_vec[i]
        #println("P, ", P)
        T = T_vec[i]
        #println("T, ", T)
        species_moles_vec = moles_mat[:, i] 
        #println("species_moles_vec, ", species_moles_vec)
        concentrations_vec = [species_moles / section_volume for species_moles in species_moles_vec]
        #println("concentrations_vec, ", concentrations_vec)
        
        common = calculate_rates_and_properties(T, P, concentrations_vec, p) #having to put the [i] doesn't make sense

        F_flux = fill(copy(species_moles_vec) .* 1.0u"1/s", N_points)
        #V_flux = fill(0.0u"m^3/s", N_points) #might be needed later
        for (j, chemical) in enumerate(p.components)
            total_concentration = sum(species_moles_vec) / section_volume

            molar_volume = Clapeyron.volume(p.clapeyron_model, P, T, species_moles_vec) / 1.0u"mol"
            
            generation_j = common.total_species_rates[chemical] * section_volume

            #concentrations_vec[i][j] * (F_flux[j-1] * molar_volume)
            if i == 1
                volumetric_flow_rate = p.initial_F_in[j] * molar_volume 
                F_flux[i][j] = concentrations_vec[j] * volumetric_flow_rate #not entirely sure how to get volumetric flow here
                println(F_flux[i][j])
                dmoles_dt[j, i] = p.initial_F_in[j] - F_flux[i][j] + generation_j
            else
                volumetric_flow_rate = F_flux[i-1][j] * molar_volume
                F_flux[i][j] = concentrations_vec[j] * volumetric_flow_rate
                dmoles_dt[j, i] = F_flux[i-1][j] - F_flux[i][j] + generation_j
            end
        end

        dP_dt[i] = 0.0u"bar/s"
        #reactor pressure will remain unchanged for now
        
        #dJacket_fluid_mass_flow_dt = (proportional_gain * error) * 
        #(integral_time * sum(p.previous_errors)) * 
        #(derivative_time * ((p.previous_errors[end] - p.previous_errors[end-1]) / (t[end] - t[end-1]))) 
        #not sure how to get the time between the previous and this run

        
        jacket_to_wall_heat_flux = (p.reactor_wall_overall_heat_transfer_coeff * section_wall_area * (T - T_jacket_vec[i]))

        if i == 1
            initial_jacket_flow_flux = (m_flow_jacket_vec[i] * p.jacket_fluid_cp * (p.jacket_fluid_inlet_temperature - T_jacket_vec[i]))
            jacket_flow_flux[i] = (m_flow_jacket_vec[i] * p.jacket_fluid_cp * (T_jacket_vec[i+1] - T_jacket_vec[i]))
            dT_jacket_dt[i] = (jacket_to_wall_heat_flux + initial_jacket_flow_flux - jacket_flow_flux[i]) / (300.0u"g" * p.jacket_fluid_cp)
        elseif i < N_points
            jacket_flow_flux[i] = (m_flow_jacket_vec[i] * p.jacket_fluid_cp * (T_jacket_vec[i+1] - T_jacket_vec[i]))
            dT_jacket_dt[i] = (jacket_to_wall_heat_flux + jacket_flow_flux[i-1] - jacket_flow_flux[i]) / (300.0u"g" * p.jacket_fluid_cp)
        else
            jacket_flow_flux[i] = jacket_flow_flux[i-1]
            #(m_flow_jacket_vec[i] * p.jacket_fluid_cp * (300.0u"K" - T_jacket_vec[i]))
            dT_jacket_dt[i] = (jacket_to_wall_heat_flux + jacket_flow_flux[i-1] - jacket_flow_flux[i]) / (300.0u"g" * p.jacket_fluid_cp)
        end
        #this will be for jacket fluid temperature 

        net_reaction_rates = common.net_reaction_rates
        delta_H_from_reactions = 0.0u"W"
        for i in eachindex(p.reactions)
            delta_H_from_reactions += p.heat_of_reactions_dict[p.reactions[i].name] * (-net_reaction_rates[p.reactions[i].name] * section_volume)
        end

        total_heat_capacity = 0.0u"J/K"
        for (j, (chemical, concentration)) in enumerate(common.concentrations_dict)
            total_heat_capacity += p.heat_capacities_dict[chemical] * species_moles_vec[j]
        end

        section_wall_area = (p.reactor_diameter * pi) * section_length

        ambient_flux = (ambient_overall_heat_transfer_coeff * section_wall_area * (p.ambient_temperature - T))

        dT_dt[i] = (-jacket_to_wall_heat_flux + ambient_flux + delta_H_from_reactions) / total_heat_capacity
    end
    return nothing
end

function pid_condition(u, t, integrator)
    return t - integrator.p.pid_state.last_time >= 0.5u"s"
end

function pid_affect!(integrator)
    p = integrator.p
    t = integrator.t

    current_temps = integrator.u.x[2]
    avg_T = sum(current_temps) / length(current_temps)

    error = avg_T - p.desired_reactor_temp
    dt = t - p.pid_state.last_time
    
    p.pid_state.integral_error += error * dt
    
    derivative = (error - p.pid_state.last_error) / dt

    
    output_flow = 0.3u"kg/s" + (p.proportional_gain * error * 1.0u"kg/(s*K)") + #Kp
                  (p.integral_time * p.pid_state.integral_error * 1.0u"kg/(s*K*s)") + #Ki
                  (p.derivative_time * derivative * 1.0u"kg*s/(s*K)" #Kd
    )
    output_flow = clamp(output_flow, 0.01u"kg/s", 5.0u"kg/s")    
    
    for i in 1:length(integrator.u.x[4])
        integrator.u.x[4][i] = output_flow
    end

    #=
    output_flow = 300.0u"K"
    output_flow -= (p.proportional_gain * error * 1.0u"1") 
    output_flow -= (p.integral_time * p.pid_state.integral_error * 1.0u"1/s")
    output_flow -= (p.derivative_time * derivative * 1.0u"s") #Kp
    output_flow = clamp(output_flow, 200.0u"K", 1000.0u"K")  

    for i in 1:length(integrator.u.x[3])
        integrator.u.x[3][i] = output_flow
    end
    =#

    p.pid_state.last_error = error
    p.pid_state.last_time = t
end
#=
DiscreteCallback(condition, affect!;
    initialize = INITIALIZE_DEFAULT,
    finalize = FINALIZE_DEFAULT,
    save_positions = (true, true),
    initializealg = nothing)=#



#The order in which reaction defs a put in the dictionary determines the order in which they will be applied

reactor_inlet_temperature = 350.13u"°C" |> u"K"

reactor_inlet_pressure = 1.0u"bar"

SMR_model = PR(["Methanol", "water", "carbon monoxide", "hydrogen", "Carbon Dioxide"])

SMR_methanol_total_moles = 0.06520168451976938u"mol"

z = [1.0, 1.3, 0.001, 0.001, 0.001]
SMR_initial_C = z .* u"mol/L"

SMR_species_total_moles = z .* SMR_methanol_total_moles

# CH3OH -> CO + 2H2
# Below is for MD
ref_T = 300.13u"°C" |> u"K"
kf_ref = 1e-1u"1/s" 
Ea_f = 90u"kJ/mol" 
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

Keq_ref = K_gibbs_free(298u"K", ref_T, 25.2u"kJ/mol", 90.7u"kJ/mol")
#formatted reference temperature, actual temperature, delta gibbs free energy of reaction, and heat of reaction

kr_ref = 0u"L^2/(mol^2*s)"
# For Ea_r, you know Ea_r = Ea_f - dH_rxn, assuming elementary.
Ea_r = 43u"kJ/mol"
kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T) 
#The Above /1.0u"bar^1" is very necessary because it causes the net_eraction_rate to not error on a dimension error (trust me!)

MD_reaction = HomogenousChemicalReaction(
    "MD",
    Dict("CH3OH" => 1),
    Dict("CO" => 1, "H2" => 2),
    ["CH3OH", "CO", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)

# CO + H2O ⇋ CO2 + H2
# Below is for WGS
ref_T = 300.13u"°C" |> u"K"
kf_ref = 0.2e-1u"L/(mol*s)" 
Ea_f = 60u"kJ/mol" #usually around 55-100 kJ/mol
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

# For reverse reaction, calculate kr from Keq. 
Keq_ref = K_gibbs_free(298u"K", ref_T, -28.63u"kJ/mol", -41.1u"kJ/mol")
#formatted reference temperature (298K), actual temperature (~350*C), delta gibbs free energy of reaction (-28.63), and heat of reaction (-41.1)

kr_ref = kf_ref / Keq_ref

Ea_r = 100.0u"kJ/mol"

kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T)

WGS_reaction = HomogenousChemicalReaction(
    "WGS",
    Dict("CO" => 1, "H2O" => 1),
    Dict("CO2" => 1, "H2" => 1),
    ["CO", "H2O", "CO2", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)
N_points = 20

reactor_desired_temperature = 350.13u"°C" |> u"K"

proportional_gain = 8.0
integral_time = 10.0
derivative_time = 0.1

initial_prev_errors = [1.0u"K"]

reactor_length = 1.0u"m"

reactor_diameter = 0.01u"m"

catalyst_packing_density = 3000.0u"kg/m^3"

ambient_overall_heat_transfer_coeff = 1.0u"W/(m^2*K)"
ambient_temperature = 300.13u"K"

reactor_wall_overall_heat_transfer_coeff = 1.0u"W/(m^2*K)"

jacket_fluid_cp = 4.184u"J/(g*K)"

heat_capacities_dict = Dict("CH3OH" => 65.2u"J/(mol*K)", 
                            "H2O" => 30.29u"J/(mol*K)", 
                            "CO" => 29.31u"J/(mol*K)", 
                            "H2" => 46.75u"J/(mol*K)", 
                            "CO2" => 36.45u"J/(mol*K)"
)

F_in = [1.0, 1.3, 0.001, 0.001, 0.001] .* 1.0u"mol/s"

jacket_fluid_inlet_temperature = 400.13u"°C" |> u"K"

pid_state = PIDState(0.0u"K", 0.0u"K*s", 0.0u"s")

SMR_parameters = PBRModel(
    [MD_reaction, WGS_reaction], #reaction dict
    SMR_model, #thermo model
    ["CH3OH", "H2O", "CO", "H2", "CO2"], #components dict
    Dict("MD" => 90.7u"kJ/mol", "WGS" => -41.1u"kJ/mol"), #heat of reactions dict 
    N_points, 
    pid_state,
    reactor_desired_temperature,
    proportional_gain, 
    integral_time, 
    derivative_time,
    initial_prev_errors, 
    reactor_length, 
    reactor_diameter,
    catalyst_packing_density,
    ambient_overall_heat_transfer_coeff, 
    ambient_temperature, 
    reactor_wall_overall_heat_transfer_coeff, 
    jacket_fluid_cp, 
    jacket_fluid_inlet_temperature,
    heat_capacities_dict,
    F_in,
)

coolant_initial_temperature = 300.13u"K"

u0_P = fill(reactor_inlet_pressure, SMR_parameters.N_points)
u0_T_rxn = fill(reactor_inlet_temperature, SMR_parameters.N_points)
u0_T_jkt = fill(coolant_initial_temperature, SMR_parameters.N_points)
u0_flow_jkt = fill(300.0u"g/s", SMR_parameters.N_points)
u0_mol_flows = repeat((SMR_species_total_moles .* 0.001), 1, SMR_parameters.N_points)

u0_mol_flows[1, 1] += SMR_species_total_moles[1]
u0_mol_flows[2, 1] += SMR_species_total_moles[2]
u0_mol_flows[3, 1] += SMR_species_total_moles[3]
u0_mol_flows[4, 1] += SMR_species_total_moles[4]
u0_mol_flows[5, 1] += SMR_species_total_moles[5]

u0 = ArrayPartition(u0_P, u0_T_rxn, u0_T_jkt, u0_flow_jkt, u0_mol_flows)

du0 = copy(u0) * 0.0u"1/s"

test = pbr_ode_system!(du0, u0, SMR_parameters, 100)

cb = DiscreteCallback(pid_condition, pid_affect!)

tspan = (0.0u"s", 10.0u"s")

prob = ODEProblem(pbr_ode_system!, u0, tspan, SMR_parameters)

@time sol = DifferentialEquations.solve(prob, Tsit5(), callback=cb)
#Rodas4(autodiff=false) is probably better but creating a sparsity pattern seems to be impossible with unitful

t_grid = sol.t
z_grid = 1:N_points

time_idx = Observable(1)

N = SMR_parameters.N_points

pressure_data = zeros(typeof(1.0u"bar"), length(sol.t), z_grid) 
temp_data = zeros(typeof(1.0u"K"), length(sol.t), z_grid) 
jacket_temp_data = zeros(typeof(1.0u"K"), length(sol.t), z_grid) 
jacket_mass_flow_data = zeros(typeof(1.0u"kg/s"), length(sol.t), z_grid) 
F_CH3OH_data = zeros(typeof(1.0u"mol"), length(sol.t), z_grid)
F_H2_data = zeros(typeof(1.0u"mol"), length(sol.t), z_grid)
F_H2O_data = zeros(typeof(1.0u"mol"), length(sol.t), z_grid)
F_CO_data = zeros(typeof(1.0u"mol"), length(sol.t), z_grid)
F_CO2_data = zeros(typeof(1.0u"mol"), length(sol.t), z_grid)

sol.u[1].x[1]

for i in eachindex(sol.t)
    pressure_data[i, :] = sol.u[i].x[1][z_grid] #pressure
    temp_data[i, :] = sol.u[i].x[2][z_grid]
    jacket_temp_data[i, :] = sol.u[i].x[3][z_grid]
    jacket_mass_flow_data[i, :] = sol.u[i].x[4][z_grid]
    F_CH3OH_data[i, :] = sol.u[i].x[5][1, :][z_grid]
    F_H2O_data[i, :] = sol.u[i].x[5][2, :][z_grid]
    F_CO_data[i, :] = sol.u[i].x[5][3, :][z_grid]
    F_H2_data[i, :] = sol.u[i].x[5][4, :][z_grid]
    F_CO2_data[i, :] = sol.u[i].x[5][5, :][z_grid]
end

pressure_data
temp_data
jacket_temp_data
jacket_mass_flow_data
F_CH3OH_data
F_H2O_data
F_H2_data
F_CO_data
F_CO2_data

# Determine the time range for the animation
desired_max_time = sol.t[end]
max_time_idx = argmin(abs.(sol.t .- desired_max_time))
time_indices_to_animate = 1:max_time_idx

fig = Figure(size = (800, 600))

title_str = @lift("Reactor Profile at t = $(round(typeof(1.0u"s"), sol.t[$time_idx], digits=1))")
ax = Axis(fig[1, 1],
          xlabel = "Reactor Length (m)",
          ylabel = "Molar Flow (mol/s)",
          title = title_str)

# 2. "Lift" the data for each species to be plotted.
# The `$` tells @lift to use the *value* of the Observable.
# When `time_idx` changes, these y-values will automatically update.
y_CH3OH = @lift(F_CH3OH_data[$time_idx, z_grid])
y_H2O = @lift(F_H2O_data[$time_idx, z_grid])
y_H2 = @lift(F_H2_data[$time_idx, z_grid])
y_CO = @lift(F_CO_data[$time_idx, z_grid])
y_CO2 = @lift(F_CO2_data[$time_idx, z_grid])
temp_reactor = @lift(temp_data[$time_idx, z_grid])


# 3. Plot the lifted (reactive) data
# The x-data (z_grid) is static, but the y-data is an Observable.
lines!(ax, z_grid, y_CH3OH, label = "CH3OH", linewidth=2)
lines!(ax, z_grid, y_H2, label = "H2", linewidth=2)
lines!(ax, z_grid, y_H2O, label = "H2O", linewidth=2)
lines!(ax, z_grid, y_CO, label = "CO", linewidth=2)
lines!(ax, z_grid, y_CO2, label = "CO2", linewidth=2)
ax_temp = Axis(fig[1, 1],
               ylabel = "Temperature (K)",
               yaxisposition = :right
)
hidexdecorations!(ax_temp)
lines!(ax_temp, temp_reactor, label = "Temperature", linewidth=2)

# Add a legend
axislegend(ax)

# Optional: Fix the y-axis limits to prevent them from jumping around during the animation
# Find the global min/max over the animated time range and add some padding.
flow_data_subset = vcat(
    F_CH3OH_data[time_indices_to_animate, z_grid],
    F_H2_data[time_indices_to_animate, z_grid],
    F_H2O_data[time_indices_to_animate, z_grid],
    F_CO_data[time_indices_to_animate, z_grid],
    F_CO2_data[time_indices_to_animate, z_grid]
)
flow_ymin = minimum(flow_data_subset)
flow_ymax = maximum(flow_data_subset)
GLMakie.ylims!(ax, flow_ymin - 0.1 * abs(flow_ymin), flow_ymax + 0.1 * abs(flow_ymax))

# For temperature:
temp_ymin = minimum(temp_data[time_indices_to_animate, :])
temp_ymax = maximum(temp_data[time_indices_to_animate, :])
GLMakie.ylims!(ax_temp, temp_ymin - 0.1 * abs(temp_ymin), temp_ymax + 0.1 * abs(temp_ymax))

# Display the initial figure
display(fig)

for i in time_indices_to_animate
    time_idx[] = i
    sleep(0.01)
end

record_sol = false
if record_sol
    framerate = 30
    rm("C://Users//wille//Desktop//Legacy-Projects-And-Code-Will-Martin//Julia Projects//interactive_PBR.mp4", force = true)
    output_file = "C://Users//wille//Desktop//Legacy-Projects-And-Code-Will-Martin//Julia Projects//interactive_PBR.mp4"
    record(fig, output_file, time_indices_to_animate; framerate = framerate) do i
        # This block is executed for each frame.
        # All we need to do is update the time index.
        time_idx[] = i
    end
end