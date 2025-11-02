using Revise

using Unitful, DifferentialEquations
using NLsolve
using Clapeyron
using Plots
using Symbolics
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

struct ChemicalReaction
    name::String
    # Use Dicts for reactants/products with species name => stoichiometric coefficient
    # For rate law, these are the powers in elementary reactions.
    # For mass balance, they are the coefficients.
    reactants::Dict{String, Int}
    products::Dict{String, Int}
    #To maintain order throughout the process
    all_chemicals::Vector{String}
    kf_A::typeof(1.0u"L/(mol*s)") # Pre-exponential factor for forward reaction
    kf_Ea::typeof(1.0u"J/mol")   # Activation energy for forward reaction
    kr_A::typeof(1.0u"L/(mol*s)") # Pre-exponential factor for reverse reaction
    kr_Ea::typeof(1.0u"J/mol")   # Activation energy for reverse reaction
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

"""Example for Steam Propane Reforming
SPR_reaction = ChemicalReaction(
    "SPR",
    Dict("C3H8" => 1, "H2O" => 3),
    Dict("CO" => 3, "H2" => 7),
    ["C3H8", "H2O", "CO", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)
end


"""

function overall_rate(reaction, model, temperature, pressure, concentrations_vec)
    total_moles = sum(concentrations_vec)
    mole_fractions = ustrip.(concentrations_vec ./ total_moles)

    activity_coeffs = activity_coefficient(model, ustrip(uconvert(u"kPa", pressure)), ustrip(uconvert(u"K", temperature)), mole_fractions)
    activities = activity_coeffs .* concentrations_vec
    
    reactants_activities = activities[1:length(reaction.reactants)]
    products_activities = activities[length(reaction.products)+1:end]

    forward_k = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
    reverse_k = arrenhius_equation_rate_constant(reaction.kr_A, reaction.kr_Ea, temperature)
    
    if typeof(reactants_activities) == Unitful.Quantity
        store_unit = unit(reactants_activities)
        reactants_activities = ustrip.(reactants_activities)
        products_activities = ustrip.(products_activities)
        return store_unit * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
    else
        return forward_k * prod(reactants_activities) - reverse_k * prod(products_activities)
    end
end

function rate_law_for_chemicals(reaction, model, temperature, pressure, concentrations_vec)
    total_moles = sum(concentrations_vec)
    mole_fractions = ustrip.(concentrations_vec ./ total_moles)

    activity_coeffs = activity_coefficient(model, ustrip(uconvert(u"kPa", pressure)), ustrip(uconvert(u"K", temperature)), mole_fractions)
    activities = activity_coeffs .* concentrations_vec

    reactants_activities = activities[1:length(reaction.reactants)]
    products_activities = activities[length(reaction.products)+1:end]

    forward_k = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
    reverse_k = arrenhius_equation_rate_constant(reaction.kr_A, reaction.kr_Ea, temperature)
    
    stoichs_and_signs = Dict()
    for chemical in reaction.all_chemicals
        if haskey(reaction.reactants, chemical)
            stoichs_and_signs[chemical] = -1 * reaction.reactants[chemical]
        else
            stoichs_and_signs[chemical] = 1 * reaction.products[chemical]
        end
    end

    rates = Dict()
    if typeof(reactants_activities) == Unitful.Quantity
        store_unit = unit(reactants_activities)
        reactants_activities = ustrip.(reactants_activities)
        products_activities = ustrip.(products_activities)

        for chemical in reaction.all_chemicals
            if haskey(reaction.reactants, chemical)
                rates[chemical] = (store_unit * (stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))))
            else 
                rates[chemical] = (store_unit * (stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))))
            end
        end
    else
        for chemical in reaction.all_chemicals
            if haskey(reaction.reactants, chemical)
                rates[chemical] = stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
            else
                rates[chemical] = stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
            end
        end
    end
    return rates
end
function heterogeneous_rate_law_for_chemicals(reaction, model, temperature, pressure, molar_flows_vec)
    total_molar_flow = sum(molar_flows_vec)
    mole_fractions = ustrip.(molar_flows_vec ./ total_molar_flow)

    activity_coeffs = activity_coefficient(model, ustrip(uconvert(u"kPa", pressure)), ustrip(uconvert(u"K", temperature)), mole_fractions)
    activities = activity_coeffs .* mole_fractions

    reactants_activities = activities[1:length(reaction.reactants)]
    products_activities = activities[length(reaction.products)+1:end]

    forward_k = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
    reverse_k = arrenhius_equation_rate_constant(reaction.kr_A, reaction.kr_Ea, temperature)
    
    stoichs_and_signs = Dict()
    for chemical in reaction.all_chemicals
        if haskey(reaction.reactants, chemical)
            stoichs_and_signs[chemical] = -1 * reaction.reactants[chemical]
        else
            stoichs_and_signs[chemical] = 1 * reaction.products[chemical]
        end
    end

    rates = Dict()
    if typeof(reactants_activities) == Unitful.Quantity
        store_unit = unit(reactants_activities)
        reactants_activities = ustrip.(reactants_activities)
        products_activities = ustrip.(products_activities)

        for chemical in reaction.all_chemicals
            if haskey(reaction.reactants, chemical)
                rates[chemical] = (store_unit * (stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))))
            else 
                rates[chemical] = (store_unit * (stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))))
            end
        end
    else
        for chemical in reaction.all_chemicals
            if haskey(reaction.reactants, chemical)
                rates[chemical] = stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
            else
                rates[chemical] = stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
            end
        end
    end
    return rates
end

function K_eq(forward_k, reverse_k)
    return forward_k / reverse_k
end

#K1 is initial equilibrum constant, returns equilibrium constant after temperature change
function van_t_hoft(K_eq1, T1, T2, heat_of_reaction)
    R_GAS = 8.314u"J/(mol*K)"
    return K_eq1 * exp((-heat_of_reaction / R_GAS) * (1/T2 - 1/T1))
end

function Q_gibbs_free(reaction, clapeyron_model, T, p, z, heat_of_reaction)
    potentials = chemical_potential(clapeyron_model, p, T, z)
    R_gas = 8.31446261815324u"J/(K*mol)"

    delta_gibbs = 0u"kJ"
    i = 0
    for chemical in reaction.all_chemicals
        if haskey(reaction.reactants, chemical) 
            delta_gibbs += potentials[i]
        else 
            delta_gibbs -= potentials[i]
        i += 1
        end
    end
    return exp((-1*(reactant_gibbs_free_energy - product_gibbs_free_energy)) / (R_gas * T))
    #return exp((heat_of_reaction - (T * (reactant_entropy - product_entropy))) / (R_gas * T))
end  
function K_gibbs_free(T_ref, T_actual, ΔG_rxn_ref, ΔH_rxn_ref)
    T_actual = uconvert(u"K", T_actual) # Ensure T is a plain number in Kelvin

    R = 8.314e-3u"kJ/(mol*K)" # Gas constant in kJ
    
    K_ref = exp(-ΔG_rxn_ref / (R * T_ref))

    ln_K_ratio = (-ΔH_rxn_ref / R) * (1/T_actual - 1/T_ref)
    
    K_T = K_ref * exp(ln_K_ratio)
    
    return ustrip(K_T) # Return a unitless equilibrium constant
end

"""
function simple_pressure_rate_law(reaction, temperature, partial_pressures_dict)
    # Calculate the k_p at the current temperature
    kp = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
    
    rates = Dict()
    i = 1
    for chemical in reaction.all_chemicals
        if haskey(reaction.reactants, chemical)
            rates[chemical] = -1 * reaction.reactants[chemical] * partial_pressures_dict[chemical] * kp 
        else
            rates[chemical] = 1 * reaction.products[chemical] * partial_pressures_dict[chemical] * kp 
        end
        i += 1
    end
    return rates
end
"""
function reversible_pressure_rate_law(reaction, temperature, partial_pressures_dict)
    # 1. Calculate forward and reverse rate constants at the current temperature
    kf = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
    kr = arrenhius_equation_rate_constant(reaction.kr_A, reaction.kr_Ea, temperature)

    # 2. Calculate the forward reaction "propensity" based on reactant partial pressures
    forward_term = 1.0u"bar^0" # Initialize with a unitless 1
    for (reactant, nu) in reaction.reactants
        # Ensure we don't take a negative pressure to a power
        pressure = max(partial_pressures_dict[reactant], 0.0u"bar")
        println("forward pressre ", pressure, partial_pressures_dict[reactant])
        forward_term *= pressure^nu
    end

    # 3. Calculate the reverse reaction "propensity" based on product partial pressures
    reverse_term = 1.0u"bar^0" # Initialize with a unitless 1
    for (product, nu) in reaction.products
        pressure = max(partial_pressures_dict[product], 0.0u"bar")
        println("reverse pressre ", pressure, partial_pressures_dict[product])
        reverse_term *= pressure^nu
    end
    
    # 4. Calculate the net rate of reaction (mol_rxn / kg_cat / s)
    net_reaction_rate = (ustrip(kf) * ustrip(forward_term) - ustrip(kr) * ustrip(reverse_term)) * 1.0u"mol/(kg*s)"
    #0.00001u"mol/(kg*s)"
    println("name: ", reaction.name, ", kf: ", kf, ", Forward Term: ", forward_term, ", kr: ", kr, ", Reverse Term: ", reverse_term, ", Net Reaction Rate: ", ustrip(kf) - ustrip(kr))
    #(ustrip(kf) - ustrip(kf)) * 1.0u"mol/(kg*s)"
    #kf * forward_term - kr * reverse_term
    
    # 5. Calculate the rate of formation/consumption for each species (mol_species / kg_cat / s)
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
        println("chemical: ", chemical, ", net_reaction_rate: ", net_reaction_rate, ", net stoichiometry: ", net_stoichiometry, ", rate for chemical: ", rates[chemical])
    end
    
    return rates
end


"""
Things to Implement
    - add WGS to main reactor
    - Implement energy balance 
        - dT/dW = (U*a*(T_ambient - T) + r_SMR*(-ΔH_SMR) + r_WGS*(-ΔH_WGS)) / (Σ(F_i * C_pi))
            - U*a is the overall heat transfer coefficient times the area per unit volume.
            - T_ambient is the temperature of the heating medium.
            - r_SMR and r_WGS are the net rates of reaction (mol/kg/s).
            - ΔH is the heat of reaction (make sure the sign is correct!). SMR is endothermic, so (-ΔH_SMR) will be a 
            negative term (cooling). WGS is exothermic.
            - F_i are the molar flows and C_pi are the molar heat capacities of each component.
    - Implement pressure drop
        - Add pressure state to reactor
        - Add the Ergun equation as an ODE for dP/dW. You will need catalyst and reactor properties like catalyst particle 
        diameter, bed void fraction, and reactor cross-sectional area.
"""

reactor_temperature = 550.13u"°C" |> u"K"

SMR_methanol_molar_flow = 0.06520168451976938u"mol/s"

SMR_model = PR(["Methanol", "water", "carbon monoxide", "hydrogen", "Carbon Dioxide"])

z = [1.0, 3.0, 0.0001, 0.0001, 0.0001]
#z = [1, 1.3, 0.1, 0.1, 0.001]
SMR_initial_C = z .* u"mol/L"

# CH3OH + H2O ⇋ CO + 3H2
#K_gibbs_free(temporary_reaction, SMR_model, reactor_temperature, 1u"atm", z, 49.4u"kJ/mol")
# Below is for SMR
ref_T = reactor_temperature
kf_ref = 1e-1u"mol/(kg*s*bar^2)" #around 1e-5 to 1e-4
Ea_f = 100u"kJ/mol" #usually around 55-100 kJ/mol
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

# For reverse reaction, calculate kr from Keq. 
Keq_ref = K_gibbs_free(298u"K", reactor_temperature, 33.2u"kJ/mol", 49.4u"kJ/mol")
#8587
#K_gibbs_free(temporary_reaction, SMR_model, reactor_temperature, 1u"atm", z, 49.4u"kJ/mol")
kr_ref = kf_ref / Keq_ref
# For Ea_r, you know Ea_r = Ea_f - dH_rxn, assuming elementary.
# dH_rxn = +33.8 kJ/mol (from our previous discussion)
Ea_r = 50u"kJ/mol"
#Ea_f - 500u"kJ/mol" # Be careful with units here! Ea_f and dH_rxn need to be in same units.
kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T) 
#/ 1.0u"bar^2"

SMR_reaction = HeterogeneousChemicalReaction(
    "SMR",
    Dict("CH3OH" => 1, "H2O" => 1),
    Dict("CO" => 1, "H2" => 3),
    ["CH3OH", "H2O", "CO", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)

#rate_law_for_chemicals(SMR_reaction, SMR_model, reactor_temperature, 100000u"Pa", SMR_initial_C)


# CO + H2O ⇋ CO2 + H2
#K_gibbs_free(temporary_reaction, SMR_model, reactor_temperature, 1u"atm", z, 49.4u"kJ/mol")
# Below is for WGS
ref_T = reactor_temperature
kf_ref = 1e-2u"mol/(kg*s*bar^2)" #around 1e-6 to 1e-5
Ea_f = 120u"kJ/mol" #usually around 55-100 kJ/mol
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

# For reverse reaction, calculate kr from Keq. 
Keq_ref = K_gibbs_free(298u"K", reactor_temperature, -28.63u"kJ/mol", -41.1u"kJ/mol")
#8587
#K_gibbs_free(temporary_reaction, SMR_model, reactor_temperature, 1u"atm", z, 49.4u"kJ/mol")
kr_ref = kf_ref / Keq_ref
# For Ea_r, you know Ea_r = Ea_f - dH_rxn, assuming elementary.
# dH_rxn = +33.8 kJ/mol (from our previous discussion)
Ea_r = 70.0u"kJ/mol"
#Ea_f - 500u"kJ/mol" # Be careful with units here! Ea_f and dH_rxn need to be in same units.
kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T)
#455.90239536593u"mol/(kg*Pa*s)"#
WGS_reaction = HeterogeneousChemicalReaction(
    "WGS",
    Dict("CO" => 1, "H2O" => 1),
    Dict("CO2" => 1, "H2" => 1),
    ["CO", "H2O", "CO2", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)

#overall_rate(SMR_reaction, SMR_model, reactor_temperature, 100000u"Pa", SMR_initial_C)

molar_flow_correction_factor = SMR_methanol_molar_flow / SMR_initial_C[1]

SMR_molar_flows = SMR_initial_C .* molar_flow_correction_factor

total_molar_flow = sum(SMR_molar_flows)


function pbr_ode_system!(dF_dW, molar_flows, p, catalyst_weight)
    reactions = p.reaction_defs
    model = p.clapeyron_model 
    temperature = p.temperature
    pressure = p.pressure 
    components = p.components
    #input_molar_flow = p.molar_flow

    total_molar_flow = sum(molar_flows)
    println("Components: ", components, " Molar Flows: ", molar_flows)
    molar_flows_dict = Dict(zip(components, molar_flows))

    partial_pressures_dict = Dict()
    for (chemical, molar_flow) in molar_flows_dict
        partial_pressures_dict[chemical] = (molar_flow / total_molar_flow) * pressure
    end

    total_species_rates = Dict(c => 0.0u"mol/s/kg" for c in components)
    for i in eachindex(reactions)
        current_partial_pressure_dict = Dict()
        for chemical in reactions[i].all_chemicals
            current_partial_pressure_dict[chemical] = (molar_flows_dict[chemical] / total_molar_flow) * pressure
        end

        current_chemical_species_rates = reversible_pressure_rate_law(reactions[i], temperature, current_partial_pressure_dict)
        println(reactions[i].name, reactions[i].all_chemicals, current_chemical_species_rates)
        for (chemical, rate) in current_chemical_species_rates
            total_species_rates[chemical] += current_chemical_species_rates[chemical]
        end
    end

    #dF_dW *= 0.0u"mol/s/kg"
    #i is for each reaction, j is for each chemical
    
    for (i, chemical) in enumerate(components)
        dF_dW[i] = total_species_rates[chemical]
    end

    println(catalyst_weight)
    """
    for i in range(1, length(reactions))
        j = 1
        for chemical in reactions[i].all_chemicals
            dF_dW[j] = rates_dict[chemical]
        j += 1
        end
    end
    """
    return nothing
end

#The order in which reaction defs are put in the dictionary determines the order in which they will be applied
SMR_parameters = (
    reaction_defs = [SMR_reaction, WGS_reaction], 
    clapeyron_model = SMR_model, 
    temperature = reactor_temperature, 
    components = ["CH3OH", "H2O", "CO", "H2", "CO2"],
    pressure = 1.0u"bar",
)

Wspan = (0.0u"kg", 200.0u"kg")

u0 = SMR_molar_flows 

pbr_ode_problem = ODEProblem(pbr_ode_system!, u0, Wspan, SMR_parameters)

pbr_sol = DifferentialEquations.solve(pbr_ode_problem, Tsit5())

mole_flow_Methanol_values = [u[1] for u in pbr_sol.u]
mole_flow_H2O_values = [u[2] for u in pbr_sol.u]
mole_flow_CO_values = [u[3] for u in pbr_sol.u]
mole_flow_H2_values = [u[4] for u in pbr_sol.u]
mole_flow_CO2_values = [u[5] for u in pbr_sol.u]

catalyst_weights = pbr_sol.t

packing_density = 0.5u"cm^3/cm^3"

catalyst_density = 8.96u"g/cm^3"

approximate_reactor_volues = catalyst_weights ./ (catalyst_density * packing_density)
approximate_reactor_volues = uconvert.(u"L", approximate_reactor_volues)

plt = plot(catalyst_weights, mole_flow_Methanol_values,
    xlabel="Catalyst Weight",
    ylabel="Molar Flow Through Reactor",
    title="PBR Species Molar Flow Through Reactor and \nHeat Input vs. Catalyst Weight",
    legend=:best, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    label="CH3OH",
    framestyle=:box)
    
# 2. Add the other species to the existing plot using plot!()
plot!(catalyst_weights, mole_flow_H2O_values, linewidth=2, label="H2O")
plot!(catalyst_weights, mole_flow_CO_values, linewidth=2, label="CO")
plot!(catalyst_weights, mole_flow_H2_values, linewidth=2, label="H2")
plot!(catalyst_weights, mole_flow_CO2_values, linewidth=2, label="CO2")


wattages_required = []
for i in eachindex(mole_flow_Methanol_values)
    push!(wattages_required, uconvert(u"W", (mole_flow_Methanol_values[1] - mole_flow_Methanol_values[i]) * 49.4u"kJ/mol"))
    wattages_required[i] += uconvert(u"W", (mole_flow_CO2_values[1] - mole_flow_CO2_values[i]) * -41.1u"kJ/mol")
end
wattages_required = ustrip.(wattages_required)

plt = plot!(twinx(plt), catalyst_weights, wattages_required,
    ylabel="Heat Input (W)",
    legend=:bottomright, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    color=:red,
    label="Heat Input (W)",
    framestyle=:box)

display(plt)
    
plt2 = plot(catalyst_weights, mole_flow_Methanol_values,
    xlabel="Catalyst Weight",
    ylabel="Molar Flow Through Reactor",
    title="PBR Species Molar Flow Through Reactor and \nHeat Input vs. Catalyst Weight",
    legend=:false, 
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    framestyle=:box)
    
# 2. Add the other species to the existing plot using plot!()
plot!(catalyst_weights, mole_flow_H2O_values, linewidth=2)
plot!(catalyst_weights, mole_flow_CO_values, linewidth=2)
plot!(catalyst_weights, mole_flow_H2_values, linewidth=2)
plot!(catalyst_weights, mole_flow_CO2_values, linewidth=2)


wattages_required = []
for i in eachindex(mole_flow_Methanol_values)
    push!(wattages_required, uconvert(u"W", (mole_flow_Methanol_values[1] - mole_flow_Methanol_values[i]) * 49.4u"kJ/mol"))
    wattages_required[i] += uconvert(u"W", (mole_flow_CO2_values[1] - mole_flow_CO2_values[i]) * -41.1u"kJ/mol")
end

wattages_required = ustrip.(wattages_required)

plt2 = plot!(twinx(plt2), catalyst_weights, wattages_required,
    ylabel="Heat Input (W)",
    legend=:false, 
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    color=:red,
    framestyle=:box)

display(plt2)











"""
function cstr_nl_system!(F, outlet_molar_flows, p)
    reaction = p.reaction_def 
    model = p.clapeyron_model 
    temperature = p.temperature
    pressure= p.pressure 
    inlet_molar_flows = p.inlet_molar_flows 
    catalyst_weight = p.catalyst_weight 

    total_molar_flow = sum(inlet_molar_flows)

    partial_pressures = []
    for molar_flow in inlet_molar_flows
        push!(partial_pressures, ((molar_flow / total_molar_flow) * pressure))
    end

    rates_dict = simple_pressure_rate_law(reaction, temperature, partial_pressures)
    i = 1
    for chemical in reaction.all_chemicals
        molar_flow_of_chem_in = inlet_molar_flows[i]
        molar_flow_of_chem_out = outlet_molar_flows[i]
        F[i] = (ustrip(molar_flow_of_chem_in) - ustrip(molar_flow_of_chem_out)) + (ustrip(rates_dict[chemical]) * ustrip(catalyst_weight))
        i += 1
    end
    return nothing
end

cstr_parameters = (
    reaction_def = SMR_reaction,
    clapeyron_model = SMR_model, 
    temperature = reactor_temperature,
    pressure = 1.0u"atm",
    inlet_molar_flows = (SMR_initial_C .* SMR_methanol_molar_flow),
    catalyst_weight = 5.0u"kg",
)


catalyst_weights = 1.0u"kg" : 50u"kg" : 2000.0u"kg" #Formatted start, step, stop
# Adjust the range and step size as appropriate for your expected conversions.

# --- Initialize arrays to store results for plotting ---
catalyst_weights_for_plot = Float64[] # To store reactor volumes (stripped for plotting)
molar_flows_for_plots = [] # To store corresponding conversions (as fractions)

println("Simulating CSTR over a range of volumes...")
for M_catalyst in catalyst_weights
    # Update the reactor_volume in the parameters for the current iteration
    current_cstr_parameters = (
        reaction_def = SMR_reaction,
        clapeyron_model = SMR_model, # Your Clapeyron model instance
        temperature = reactor_temperature, # Isothermal assumption for now
        pressure = 1.0u"atm", # Constant pressure
        inlet_molar_flows = (SMR_initial_C .* SMR_methanol_molar_flow),
        catalyst_weight = M_catalyst,
    )

    initial_guess_molar_flows_stipped = ustrip(SMR_molar_flows)

    # Call nlsolve with the updated parameters
    cstr_sol = nlsolve((F, x) -> cstr_nl_system!(F, x, current_cstr_parameters), initial_guess_molar_flows_stipped)

    if converged(cstr_sol)
        # Extract outlet concentrations
        outlet_molar_flows_stripped = cstr_sol.zero

        #Keep the Unitful.unit
        outlet_molar_flows_with_units = outlet_molar_flows_stripped .* Unitful.unit(SMR_molar_flows[1]) # Re-apply units

        # Store the results for plotting
        push!(catalyst_weights_for_plot, ustrip(M_catalyst))
        push!(molar_flows_for_plots, ustrip(outlet_molar_flows_with_units[1]))# Conversion is dimensionless
    else
        # Handle non-convergence (e.g., print a warning, skip this point)
        println("Warning: Solver did not converge for V_reactor")
        # You might want to break here or try a different initial guess/solver options
    end
    # Optional: Update initial_guess_C_out_stripped for the next iteration
    # This can help convergence for the next point if volumes are increasing gradually.
    initial_guess_molar_flows_stipped = cstr_sol.zero
end

catalyst_weights_for_plot
molar_flows_for_plots = molar_flows_for_plots .* 1.0u"mol/s"

plt = plot(catalyst_weights_for_plot, molar_flows_for_plots,
           xlabel="Reactor Volume (L)",
           ylabel="Methanol Molar Flow",
           title="CSTR Molar Flow vs. Volume and Heat Required",
           legend=:false, # No legend needed for a single line
           linewidth=2,
           marker=:circle, # Add markers for data points
           markersize=3,
           grid=true,
           label="Methanol Molar Flow",
           framestyle=:box)

wattages_for_plot = molar_flows_for_plots * 49.4u"kJ/mol"

wattages_for_plot = uconvert.(u"W", wattages_for_plot)

plt = plot!(twinx(plt), catalyst_weights_for_plot, wattages_for_plot,
           ylabel="Heat Input",
           legend=:bottomright, # No legend needed for a single line
           linewidth=2,
           marker=:circle, # Add markers for data points
           markersize=3,
           grid=true,
           label="Methanol Molar Flow and Heat Input",
           framestyle=:box)


#Below is for WGS
ref_T = reactor_temperature
kf_ref = 8.0e-5u"L/(mol*s)"
Ea_f = 60u"kJ/mol"
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

# For reverse reaction, calculate kr from Keq. Let's assume Keq = 4 at 300K
Keq_ref = 4.0u"NoUnits" # Keq is dimensionless
kr_ref = kf_ref / Keq_ref
# For Ea_r, you know Ea_r = Ea_f - dH_rxn, assuming elementary.
# dH_rxn = +33.8 kJ/mol (from our previous discussion)
Ea_r = Ea_f - (-42.1)u"kJ/mol" # Be careful with units here! Ea_f and dH_rxn need to be in same units.
kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T)

# CO + H2O ⇋ CO2 + H2
WGS_reaction = ChemicalReaction(
    "WGS",
    Dict("CO" => 1, "H2O" => 1),
    Dict("CO2" => 1, "H2" => 1),
    ["CO", "H2O", "CO2", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)

WGS_model = PR(["Carbon Monoxide", "Water", "Carbon Dioxide", "Hydrogen"])

#carbon monoxide, water, carbon dioxide, hydrogen
WGS_initial_C = [pbr_sol.u[end][3], pbr_sol.u[end][2], 0.0001u"mol/L", pbr_sol.u[end][4]]  #[3, 1, 0.1, 7]

#rate_law_for_chemicals(WGS_reaction, WGS_model, reactor_temperature, 100000u"Pa", WGS_initial_C)


@variables N_H2, N_CO, N_CO2, N_CH3OH, N_H2O_prod
@variables S_to_C_ratio, O_to_C_ratio, Carbon_balance
@variables K_WGS

mole_balances = [
    3 * N_H2 ~ 1 * N_CO,
    Carbon_balance ~ N_CO + N_CO2,
    4 + 2 * S_to_C_ratio ~ 2 * N_H2 + 2 * N_H2O_prod, #4 * N_CH3OH + 2 * S_to_C_ratio ~ N_CO + 2*N_CO2 + N_H2O_prod,
    1 + 1 * O_to_C_ratio + S_to_C_ratio ~ N_CO + 2 * N_CO2 + N_H2O_prod,
    K_WGS ~ (N_CO2 * N_H2) / (N_CO * N_H2O_prod)
]


known_values = Dict(
    K_WGS          => 0.001,
    S_to_C_ratio   => 1.3
)

# 4. Variables to solve for
#    We want to find N_H2, N_CO, N_CO2, N_H2O_prod, and O_to_C_ratio
variables_to_solve = [N_H2, N_CO, N_CO2, N_H2O_prod, O_to_C_ratio]

# Step 1: Substitute known values
eqs_subs = Symbolics.substitute.(mole_balances, (known_values,))

# Step 2: Solve N_H2 and N_CO2 in terms of N_CO
sol_N_H2 = symbolic_linear_solve(eqs_subs[1], N_H2) # N_H2 ~ 2*N_CO
sol_N_CO2 = symbolic_linear_solve(eqs_subs[2], N_CO2) # N_CO2 ~ 1 - N_CO

# Substitute N_H2 and N_CO2 into the WGS equation to get N_H2O_prod in terms of N_CO
wgs_eq_simplified = Symbolics.substitute(eqs_subs[5], Dict(N_H2 => sol_N_H2[], N_CO2 => sol_N_CO2[]))
# K_WGS ~ ((1 - N_CO) * (2 * N_CO)) / (N_CO * N_H2O_prod)
# K_WGS ~ 2 * (1 - N_CO) / N_H2O_prod
# N_H2O_prod ~ 2 * (1 - N_CO) / K_WGS
sol_N_H2O_prod = 2 * (1 - N_CO) / known_values[K_WGS]

# Substitute all into the H Balance to solve for N_CO
h_balance_simplified = Symbolics.substitute(eqs_subs[3], Dict(N_H2 => sol_N_H2[], N_CO2 => sol_N_CO2[], N_H2O_prod => sol_N_H2O_prod[]))
# 4 + 2 * 0.8 ~ 2 * (2 * N_CO) + 2 * (2 * (1 - N_CO) / 0.364)
# 5.6 ~ 4 * N_CO + 10.989 * (1 - N_CO)
sol_N_CO = symbolic_linear_solve(h_balance_simplified, N_CO)
# Now we have N_CO. Let's evaluate it numerically
N_CO_val = Symbolics.value(sol_N_CO[])

# Calculate the rest of the values numerically
N_H2_val = 2 * N_CO_val
N_CO2_val = 1 - N_CO_val
N_H2O_prod_val = 2 * (1 - N_CO_val) / Symbolics.value(known_values[K_WGS])

# Finally, substitute all knowns into the O Balance to solve for O_to_C_ratio
o_balance_simplified = Symbolics.substitute(eqs_subs[4], Dict(N_CO => N_CO_val, N_CO2 => N_CO2_val, N_H2O_prod => N_H2O_prod_val))
sol_O_to_C_ratio = symbolic_linear_solve(o_balance_simplified, O_to_C_ratio)
O_to_C_ratio_val = Symbolics.value(sol_O_to_C_ratio[])

# Print results
println("--- Calculated Molar Flows (per 1 mole CH3OH feed) ---")
println("N_CO: ", N_CO_val); N_CO_val
println("N_H2: ", N_H2_val); N_H2_val
println("N_CO2: ", N_CO2_val); N_CO2_val
println("N_H2O_prod: ", N_H2O_prod_val); N_H2O_prod_val
println("-----------------------------------------------------")
println("S/C Ratio (input): ", Symbolics.value(known_values[S_to_C_ratio])); Symbolics.value(known_values[S_to_C_ratio])
println("O/C Ratio (output): ", O_to_C_ratio_val); O_to_C_ratio_val
"""




