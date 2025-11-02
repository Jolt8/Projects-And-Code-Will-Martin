using Revise

using Unitful, DifferentialEquations
using NLsolve
using Clapeyron
using Plots
using Symbolics

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
    kf_A::typeof(1.0u"mol/(kg*s*Pa)") # Pre-exponential factor for forward reaction
    kf_Ea::typeof(1.0u"J/mol")   # Activation energy for forward reaction
    kr_A::typeof(1.0u"mol/(kg*s*Pa)") # Pre-exponential factor for reverse reaction
    kr_Ea::typeof(1.0u"J/mol")   # Activation energy for reverse reaction
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
 # module Chemeqs


function simple_pressure_rate_law(reaction, temperature, partial_pressures_dict)
    # Calculate the k_p at the current temperature
    kp = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
    
    rates = Dict()
    i = 1
    for chemical in reaction.all_chemicals
        if haskey(reaction.reactants, chemical)
            rates[chemical] = -1 * reaction.reactants[chemical] * partial_pressures_dict[i] * kp 
        else
            rates[chemical] = 1 * reaction.products[chemical] * partial_pressures_dict[i] * kp 
        end
        i += 1
    end
    return rates
end







reactor_temperature = 350.13u"°C" |> u"K"

SMR_methanol_molar_flow = 0.06520168451976938u"mol/s"

temporary_reaction = ChemicalReaction(
    "SMR",
    Dict("CH3OH" => 1, "H2O" => 1),
    Dict("CO" => 1, "H2" => 3),
    ["CH3OH", "H2O", "CO", "H2"],
    1.0u"L/(mol*s)", 1.0u"kJ/mol",
    1.0u"L/(mol*s)", 1.0u"kJ/mol"
)

SMR_model = PR(["Methanol", "water", "carbon monoxide", "hydrogen"])

z = [1, 1.3, 0.001, 0.001]
SMR_initial_C = z .* u"mol/L"

# CH3OH + H2O ⇋ CO + 3H2
#K_gibbs_free(temporary_reaction, SMR_model, reactor_temperature, 1u"atm", z, 49.4u"kJ/mol")
# Below is for SMR
ref_T = reactor_temperature
kf_ref = 0.01u"mol/(kg*s*Pa)"
Ea_f = 70u"kJ/mol" #usually around 55-100 kJ/mol
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

# For reverse reaction, calculate kr from Keq. 
Keq_ref = K_gibbs_free(298u"K", reactor_temperature, 3.3u"kJ/mol", 49.4u"kJ/mol")
#8587
#K_gibbs_free(temporary_reaction, SMR_model, reactor_temperature, 1u"atm", z, 49.4u"kJ/mol")
kr_ref = kf_ref / Keq_ref
# For Ea_r, you know Ea_r = Ea_f - dH_rxn, assuming elementary.
# dH_rxn = +33.8 kJ/mol (from our previous discussion)
Ea_r = 20.6u"kJ/mol"
#Ea_f - 500u"kJ/mol" # Be careful with units here! Ea_f and dH_rxn need to be in same units.
kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T)

SMR_reaction = HeterogeneousChemicalReaction(
    "SMR",
    Dict("CH3OH" => 1, "H2O" => 1),
    Dict("CO" => 1, "H2" => 3),
    ["CH3OH", "H2O", "CO", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)


#overall_rate(SMR_reaction, SMR_model, reactor_temperature, 100000u"Pa", SMR_initial_C)

rate_law_for_chemicals(SMR_reaction, SMR_model, reactor_temperature, 100000u"Pa", SMR_initial_C)

molar_flow_correction_factor = SMR_methanol_molar_flow / SMR_initial_C[1]

SMR_molar_flows = SMR_initial_C .* molar_flow_correction_factor

total_molar_flow = sum(SMR_molar_flows)


function pbr_ode_system!(dC_dW, molar_flows, p, catalyst_weight)
    reaction = p.reaction_def 
    model = p.clapeyron_model 
    temperature = p.temperature
    pressure = p.pressure 
    input_molar_flow = p.molar_flow

    total_molar_flow = sum(molar_flows)

    partial_pressures = []
    for molar_flow in molar_flows
        push!(partial_pressures, ((molar_flow / total_molar_flow) * pressure))
    end
    #partial_pressures = (molar_flows ./ total_molar_flow) * pressure

    println(partial_pressures)

    rates_dict = simple_pressure_rate_law(reaction, temperature, partial_pressures)

    println(rates_dict)

    i = 1
    for chemical in reaction.all_chemicals
        dC_dW[i] = rates_dict[chemical]
       i += 1
    end
    return nothing
end

SMR_parameters = (
    reaction_def = SMR_reaction, 
    clapeyron_model = SMR_model, 
    temperature = reactor_temperature, 
    pressure = 100000u"Pa",
    molar_flow = SMR_methanol_molar_flow
)

Wspan = (0.0u"kg", 0.00025u"kg")

u0 = SMR_molar_flows

pbr_ode_problem = ODEProblem(pbr_ode_system!, u0, Wspan, SMR_parameters)

pbr_sol = DifferentialEquations.solve(pbr_ode_problem, Tsit5())

mole_flow_Methanol_values = [u[1] for u in pbr_sol.u]
mole_flow_H2O_values = [u[2] for u in pbr_sol.u]
mole_flow_CO_values = [u[3] for u in pbr_sol.u]
mole_flow_H2_values = [u[4] for u in pbr_sol.u]

catalyst_weights = pbr_sol.t

packing_density = 0.5u"cm^3/cm^3"

catalyst_density = 8.96u"g/cm^3"

approximate_reactor_volues = catalyst_weights ./ (catalyst_density * packing_density)
approximate_reactor_volues = uconvert.(u"L", approximate_reactor_volues)

plt = plot(catalyst_weights, mole_flow_Methanol_values,
    xlabel="Catalyst Weight",
    ylabel="Methanol Molar Flow Through Reactor",
    title="PBR Methnaol Molar Flow Through Reactor and \nHeat Input vs. Catalyst Weight",
    legend=:topleft, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    label="Methanol Molar Flow Through Reactor (mol/L)",
    framestyle=:box)

wattages_required = []
for i in range(1, length(mole_flow_Methanol_values))
    push!(wattages_required, uconvert(u"W", ((mole_flow_Methanol_values[1] - mole_flow_Methanol_values[i]) / mole_flow_Methanol_values[1]) * 1 * SMR_methanol_molar_flow * 49.4u"kJ/mol") )
end
wattages_required = ustrip.(wattages_required)

plt = plot!(twinx(plt), catalyst_weights, wattages_required,
    ylabel="Heat Input (W)",
    legend=:bottomleft, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    color=:red,
    label="Heat Input (W)",
    framestyle=:box)








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
