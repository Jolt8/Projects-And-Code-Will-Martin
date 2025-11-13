using Unitful, DifferentialEquations
using Clapeyron
using Plots

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
    println("name: ", reaction.name, ", kf: ", kf, ", Forward Term: ", forward_term, ", kr: ", kr, ", Reverse Term: ", reverse_term, ", Net Reaction Rate: ", (kf * forward_term) - (kr * reverse_term))
    println("")
    
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

abstract type ThermalModel end
abstract type PressureModel end

struct Isothermal <: ThermalModel
    reactor_temperature::typeof(1.0u"K")
end

struct NonIsothermal <: ThermalModel
    ambient_temperature::typeof(1.0u"K")
    overall_heat_transfer_coeff::typeof(1.0u"W/(m^2*K)")
    reactor_surface_area::typeof(1.0u"m^2/kg")
end

struct InputWattage <: ThermalModel
    input_wattage_per_kg::typeof(1.0u"W/kg")
end

struct CSTRModel{T<:ThermalModel}
    thermal_model::T
    reactor_pressure::typeof(1.0u"bar")
    reactor_volume::typeof(1.0u"l")
    v_out::typeof(1.0u"l/s")
    F_in::typeof([1.0u"mol/s"])
    reactions::Vector{HomogenousChemicalReaction}
    clapeyron_model::Any
    components::Vector{String}
    heat_of_reactions_dict::Dict{String, typeof(1.0u"kJ/mol")}
end


function calculate_rates_and_properties(T::typeof(1.0u"K"), P, concentrations_vec, p::CSTRModel)
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

        println(reactions[i].name, reactions[i].all_chemicals, species_rates)
        for (chemical, rate) in species_rates
            total_species_rates[chemical] += rate
        end
    end
    return (net_reaction_rates=net_reaction_rates, total_species_rates=total_species_rates, concentrations_dict=concentrations_dict, z_vec=z_vec)
end

function cstr_ode_system!(du, u, p::CSTRModel{Isothermal}, time)
    T = p.thermal_model.reactor_temperature
    P = p.reactor_pressure
    reactor_volume = p.reactor_volume
    v_out = p.v_out
    F_in = p.F_in
    println("Next Iter")
    println("T, ", T)
    println("P, ", P)
    concentrations_vec = [moles / reactor_volume for moles in u]

    println(concentrations_vec)
    
    common = calculate_rates_and_properties(T, P, concentrations_vec, p)

    for (i, chemical) in enumerate(p.components)
        F_out_i = concentrations_vec[i] * v_out

        generation_i = common.total_species_rates[chemical] * reactor_volume
        
        du[i] = F_in[i] - F_out_i + generation_i
    end
    return nothing
end

#The order in which reaction defs a put in the dictionary determines the order in which they will be applied

reactor_initial_temperature = 350.13u"°C" |> u"K"

reactor_pressure = 1.0u"bar"

reactor_volume = 10.0u"L"

v_out = 1.0u"L/s"

SMR_model = PR(["Methanol", "water", "carbon monoxide", "hydrogen", "Carbon Dioxide"])

SMR_methanol_total_moles = 0.06520168451976938u"mol"

z = [1.0, 1.3, 0.001, 0.001, 0.001]
SMR_initial_C = z .* u"mol/L"

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

SMR_moles = SMR_initial_C .* reactor_volume

tspan = (0.0u"s", 300u"s")

thermal_model = Isothermal(
    reactor_initial_temperature
)

F_in = [1.0, 1.3, 0.001, 0.001, 0.001] .* 1.0u"mol/s"

SMR_parameters = CSTRModel(
    thermal_model, #heat solver type
    reactor_pressure,
    reactor_volume,
    v_out,
    F_in,
    [MD_reaction], #reaction dict
    SMR_model, #clapeyron model
    ["CH3OH", "H2O", "CO", "H2", "CO2"], #components dict
    Dict("MD" => 90.7u"kJ/mol"), #heat of reactions dict 
)

if typeof(SMR_parameters.thermal_model) == Isothermal
    u0 = vcat(SMR_moles)
else
    u0 = vcat(reactor_initial_temperature, SMR_moles)
end

cstr_ode_problem = ODEProblem(cstr_ode_system!, u0, tspan, SMR_parameters)

cstr_sol = DifferentialEquations.solve(cstr_ode_problem, Tsit5())

u_start = length(u0) - length(SMR_parameters.components) + 1

if typeof(SMR_parameters.thermal_model) == Isothermal
    mole_flow_Methanol_values = [(u[1] / reactor_volume) * v_out for u in cstr_sol.u]
    mole_flow_H2O_values = [(u[2] / reactor_volume) * v_out for u in cstr_sol.u]
    mole_flow_CO_values = [(u[3] / reactor_volume) * v_out for u in cstr_sol.u]
    mole_flow_H2_values = [(u[4] / reactor_volume) * v_out for u in cstr_sol.u]
    mole_flow_CO2_values = [(u[5] / reactor_volume) * v_out for u in cstr_sol.u]
else
    temperature_over_reactor = [u[1] for u in cstr_sol.u]

    mole_flow_Methanol_values = [(u[2] / reactor_volume) * v_out for u in cstr_sol.u]
    mole_flow_H2O_values = [(u[3] / reactor_volume) * v_out for u in cstr_sol.u]
    mole_flow_CO_values = [(u[4] / reactor_volume) * v_out for u in cstr_sol.u]
    mole_flow_H2_values = [(u[5] / reactor_volume) * v_out for u in cstr_sol.u]
    mole_flow_CO2_values = [(u[6] / reactor_volume) * v_out for u in cstr_sol.u]
end

timestamps = cstr_sol.t

plt = plot(timestamps, mole_flow_Methanol_values,
    xlabel="Time (s)",
    ylabel="Molar Flow Through Reactor",
    title="CSTR Species Molar Flow Through \n Reactor vs. Time at 350°C",
    legend=:best, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    label="CH3OH",
    framestyle=:box,
    dpi=1000)
    
# 2. Add the other species to the existing plot using plot!()
plot!(timestamps, mole_flow_H2O_values, linewidth=2, marker=:circle, markersize=3, label="H2O")
plot!(timestamps, mole_flow_CO_values, linewidth=2, marker=:circle, markersize=3, label="CO")
plot!(timestamps, mole_flow_H2_values, linewidth=2, marker=:circle, markersize=3, label="H2")
plot!(timestamps, mole_flow_CO2_values, linewidth=2, marker=:circle, markersize=3, label="CO2")

plt2 = plot(timestamps, mole_flow_Methanol_values,
    xlabel="Time (s)",
    ylabel="Molar Flow Through Reactor",
    title="PBR Species Molar Flow Through Reactor \nvs. Time",
    legend=:false, 
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    framestyle=:box)
    
# 2. Add the other species to the existing plot using plot!()
plot!(timestamps, mole_flow_H2O_values, linewidth=2, marker=:circle, markersize=3)
plot!(timestamps, mole_flow_CO_values, linewidth=2, marker=:circle, markersize=3)
plot!(timestamps, mole_flow_H2_values, linewidth=2, marker=:circle, markersize=3)
plot!(timestamps, mole_flow_CO2_values, linewidth=2, marker=:circle, markersize=3)
