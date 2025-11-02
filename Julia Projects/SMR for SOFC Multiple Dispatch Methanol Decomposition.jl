using Unitful, DifferentialEquations
using NLsolve
using Clapeyron
using Plots
#using Symbolics
#using NaNMath
#include("Chemeqs/src/Chemeqs.jl")
#using .Chemeqs 
#import .Chemeqs:overall_rate, arrenhius_equation_pre_exponential_factor, rate_law_for_chemicals #For some reason this only works if I do this
isobaric_heat_capacity(PR("water"), 1u"bar", 573.28u"K", [1.0])


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

function reversible_pressure_rate_law(reaction, temperature, partial_pressures_dict)
    # 1. Calculate forward and reverse rate constants at the current temperature
    kf = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
    kr = arrenhius_equation_rate_constant(reaction.kr_A, reaction.kr_Ea, temperature)

    # 2. Calculate the forward reaction "propensity" based on reactant partial pressures
    forward_term = 1.0u"bar^0" # Initialize with a unitless 1
    for (reactant, nu) in reaction.reactants
        # Ensure we don't take a negative pressure to a power
        pressure = max(partial_pressures_dict[reactant], 0.0u"bar")
        println("reactant: ", reactant)
        println("pressure: ", pressure)
        println("actual pressure: ", partial_pressures_dict[reactant])
        println("nu: ", nu)
        forward_term *= pressure^nu
        println("forward term: ", forward_term)
        println("")
    end

    # 3. Calculate the reverse reaction "propensity" based on product partial pressures
    reverse_term = 1.0u"bar^0" # Initialize with a unitless 1
    for (product, nu) in reaction.products
        pressure = max(partial_pressures_dict[product], 0.0u"bar")
        #println("reverse pressre ", pressure, partial_pressures_dict[product])
        reverse_term *= pressure^nu
    end
    
    # 4. Calculate the net rate of reaction (mol_rxn / kg_cat / s)
    net_reaction_rate = ((kf * forward_term) - (kr * reverse_term))
    #0.00001u"mol/(kg*s)"
    println("name: ", reaction.name, ", kf: ", kf, ", Forward Term: ", forward_term, ", kr: ", kr, ", Reverse Term: ", reverse_term, ", Net Reaction Rate: ", (kf * forward_term) - (kr * reverse_term))
    println("")
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
        #println("chemical: ", chemical, ", net_reaction_rate: ", net_reaction_rate, ", net stoichiometry: ", net_stoichiometry, ", rate for chemical: ", rates[chemical])
    end
    
    return [rates, net_reaction_rate]
end

reactor_inlet_temperature = 300.13u"°C" |> u"K"

reactor_inlet_pressure = 1.0u"bar"

SMR_methanol_molar_flow = 0.06520168451976938u"mol/s"

SMR_model = PR(["Methanol", "water", "carbon monoxide", "hydrogen", "Carbon Dioxide"])

z = [1.0, 1.3, 0.001, 0.001, 0.001]
#z = [1, 1.3, 0.1, 0.1, 0.001]
SMR_initial_C = z .* u"mol/L"

# CH3OH -> CO + 2H2
# Below is for MD
ref_T = reactor_inlet_temperature
kf_ref = 1e-2u"mol/(kg*s*bar)" 
Ea_f = 90u"kJ/mol" 
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

Keq_ref = K_gibbs_free(298u"K", reactor_inlet_temperature, 33.2u"kJ/mol", 49.4u"kJ/mol")
#formatted reference temperature, actual temperature, delta gibbs free energy of reaction, and heat of reaction

kr_ref = 0u"mol/(kg*s*bar^3)"
# For Ea_r, you know Ea_r = Ea_f - dH_rxn, assuming elementary.
Ea_r = 43u"kJ/mol"
kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T) 
#The Above /1.0u"bar^1" is very necessary because it causes the net_eraction_rate to not error on a dimension error (trust me!)


MD_reaction = HeterogeneousChemicalReaction(
    "MD",
    Dict("CH3OH" => 1),
    Dict("CO" => 1, "H2" => 2),
    ["CH3OH", "CO", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)

# CO + H2O ⇋ CO2 + H2
# Below is for WGS
ref_T = reactor_inlet_temperature
kf_ref = 1e-2u"mol/(kg*s*bar^2)" #around 1e-6 to 1e-5
Ea_f = 60u"kJ/mol" #usually around 55-100 kJ/mol
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

# For reverse reaction, calculate kr from Keq. 
Keq_ref = K_gibbs_free(298u"K", reactor_inlet_temperature, -28.63u"kJ/mol", -41.1u"kJ/mol")
#formatted reference temperature (298K), actual temperature (~350*C), delta gibbs free energy of reaction (-28.63), and heat of reaction (-41.1)

kr_ref = kf_ref / Keq_ref

Ea_r = 100.0u"kJ/mol"

kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T)

WGS_reaction = HeterogeneousChemicalReaction(
    "WGS",
    Dict("CO" => 1, "H2O" => 1),
    Dict("CO2" => 1, "H2" => 1),
    ["CO", "H2O", "CO2", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)

test = uconvert(u"W/kg", (-90.7u"kJ/mol" * 0.02u"mol/(kg*s)") + (-(-41.1u"kJ/mol" * 0.00038896043u"mol/(kg*s)")))

#overall_rate(SMR_reaction, SMR_model, reactor_inlet_temperature, 100000u"Pa", SMR_initial_C)

molar_flow_correction_factor = SMR_methanol_molar_flow / SMR_initial_C[1]

SMR_molar_flows = SMR_initial_C .* molar_flow_correction_factor

total_molar_flow = sum(SMR_molar_flows)


abstract type ThermalModel end
abstract type PressureModel end

struct NonIsothermal <: ThermalModel
    ambient_temperature::typeof(1.0u"K")
    overall_heat_transfer_coeff::typeof(1.0u"W/(m^2*K)")
    reactor_surface_area::typeof(1.0u"m^2/kg")
end

struct NonIsothermalHX <: ThermalModel
    ambient_temperature::typeof(1.0u"K")
    overall_heat_transfer_coeff::typeof(1.0u"W/(m^2*K)")
    reactor_surface_area::typeof(1.0u"m^2/kg")
    input_wattage_per_kg::typeof(1.0u"W/kg")
end

struct Isothermal <: ThermalModel
    reactor_temperature::typeof(1.0u"K")
end

struct ErgunPressureDrop <: PressureModel
    catalyst_particle_diameter::typeof(1.0u"m")
    bed_void_fraction::Float64
    catalyst_density::typeof(1.0u"kg/m^3")
    reactor_cross_sectional_area::typeof(1.0u"m^2")
end

struct Isobaric <: PressureModel 
    reactor_pressure::typeof(1.0u"bar")
end 

struct PBRModel{T<:ThermalModel, P<:PressureModel}
    thermal_model::T
    pressure_model::P
    reactions::Vector{HeterogeneousChemicalReaction}
    clapeyron_model::PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}
    components::Vector{String}
    heat_of_reactions_dict::Dict{String, typeof(1.0u"kJ/mol")}
end


function calculate_rates_and_properties(T::typeof(1.0u"K"), P, molar_flows_vec, p::PBRModel)
    reactions = p.reactions
    components = p.components

    total_molar_flow = sum(molar_flows_vec)
    molar_flows_dict = Dict(zip(components, molar_flows_vec))
    z_vec = [molar_flow / total_molar_flow for molar_flow in molar_flows_vec]
    println("Molar flows dict: ", molar_flows_dict)

    partial_pressures_dict = Dict()
    for (chemical, molar_flow) in molar_flows_dict
        partial_pressures_dict[chemical] = (molar_flow / total_molar_flow) * P
    end
    
    net_reaction_rates = Dict()
    total_species_rates = Dict(c => 0.0u"mol/s/kg" for c in components)
    for i in eachindex(reactions)
        current_partial_pressures_dict = Dict()
        for chemical in reactions[i].all_chemicals
            current_partial_pressures_dict[chemical] = partial_pressures_dict[chemical] #this does not cause any issues
        end
        
        species_rates, net_rate = reversible_pressure_rate_law(reactions[i], T, current_partial_pressures_dict)

        net_reaction_rates[reactions[i].name] = net_rate

        println(reactions[i].name, reactions[i].all_chemicals, species_rates)
        for (chemical, rate) in species_rates
            total_species_rates[chemical] += rate
        end
    end
    return (net_reaction_rates=net_reaction_rates, total_species_rates=total_species_rates, molar_flows_dict=molar_flows_dict, molar_flows_vec=molar_flows_vec, z_vec=z_vec)
end

function calculate_dT_dW(u, p::PBRModel{NonIsothermal}, T::typeof(1.0u"K"), P, common)
    clapeyron_model = p.clapeyron_model
    components = p.components
    reactions = p.reactions
    heat_of_reactions_dict = p.heat_of_reactions_dict
    U = p.thermal_model.overall_heat_transfer_coeff
    a = p.thermal_model.reactor_surface_area
    ambient_temperature = p.thermal_model.ambient_temperature

    molar_flows_dict = common.molar_flows_dict
    net_reaction_rates = common.net_reaction_rates
    
    heat_exchange_with_environment = U*a*(ambient_temperature - T)

    delta_H_from_reactions = 0.0u"W/kg"
    for i in eachindex(reactions)
        println("current_heat_of_reaction: ", heat_of_reactions_dict[reactions[i].name])
        println("-net_reaction_rates: ", -net_reaction_rates[reactions[i].name])
        delta_H_from_reactions += heat_of_reactions_dict[reactions[i].name] * -net_reaction_rates[reactions[i].name]
    end

    total_heat_capacity_flow_rate = 0.0u"J/(K*s)"
    i = 1
    for component in components
        current_model = PR([clapeyron_model.components[i]])
        comp_heat_capacity = isobaric_heat_capacity(current_model, P, T, [1.0]) #Not sure if I should use isobaric or isochoric
        comp_heat_capacity /= 1.0u"mol" #for some reason clapeyron doesn't return this in moles just j/k
        total_heat_capacity_flow_rate += molar_flows_dict[component] * comp_heat_capacity 
        i += 1
    end

    
    println("heat_exchange_with_environment: ", heat_exchange_with_environment)
    println("delta_H_from_reactions: ", delta_H_from_reactions)
    println("total_heat_capacity_flow_rate: ", total_heat_capacity_flow_rate)
    dT_dW = (heat_exchange_with_environment + delta_H_from_reactions) / total_heat_capacity_flow_rate
    println("dT_dW: ", uconvert(u"K/kg", dT_dW))
    return dT_dW * 0.0001 #For the change in heat throughout the reactor
end

function calculate_dT_dW(u, p::PBRModel{NonIsothermalHX}, T::typeof(1.0u"K"), P, common)
    clapeyron_model = p.clapeyron_model
    components = p.components
    reactions = p.reactions
    heat_of_reactions_dict = p.heat_of_reactions_dict
    U = p.thermal_model.overall_heat_transfer_coeff
    a = p.thermal_model.reactor_surface_area
    ambient_temperature = p.thermal_model.ambient_temperature
    input_wattage_per_kg = p.thermal_model.input_wattage_per_kg

    molar_flows_dict = common.molar_flows_dict
    net_reaction_rates = common.net_reaction_rates
    
    heat_exchange_with_environment = U*a*(ambient_temperature - T)

    delta_H_from_reactions = 0.0u"W/kg"
    for i in eachindex(reactions)
        println("current_heat_of_reaction: ", heat_of_reactions_dict[reactions[i].name])
        println("-net_reaction_rates: ", -net_reaction_rates[reactions[i].name])
        delta_H_from_reactions += heat_of_reactions_dict[reactions[i].name] * -net_reaction_rates[reactions[i].name]
    end

    total_heat_capacity_flow_rate = 0.0u"J/(K*s)"
    i = 1
    for component in components
        current_model = PR([clapeyron_model.components[i]])
        comp_heat_capacity = isobaric_heat_capacity(current_model, P, T, [1.0]) #Not sure if I should use isobaric or isochoric
        comp_heat_capacity /= 1.0u"mol" #for some reason clapeyron doesn't return this in moles just j/k
        total_heat_capacity_flow_rate += molar_flows_dict[component] * comp_heat_capacity 
        i += 1
    end

    
    println("heat_exchange_with_environment: ", heat_exchange_with_environment)
    println("delta_H_from_reactions: ", delta_H_from_reactions)
    println("total_heat_capacity_flow_rate: ", total_heat_capacity_flow_rate)
    dT_dW = (heat_exchange_with_environment + delta_H_from_reactions + input_wattage_per_kg) / total_heat_capacity_flow_rate
    println("dT_dW: ", uconvert(u"K/kg", dT_dW))
    return dT_dW #For the change in heat throughout the reactor
end

function calculate_dP_dW(u, p::PBRModel{Isothermal, ErgunPressureDrop}, T::typeof(1.0u"K"), P, common)
    clapeyron_model = p.clapeyron_model
    components = p.components
    reactor_cross_sectional_area = p.pressure_model.reactor_cross_sectional_area
    bed_void_fraction = p.pressure_model.bed_void_fraction
    catalyst_particle_diameter = p.pressure_model.catalyst_particle_diameter
    catalyst_density = p.pressure_model.catalyst_density 

    molar_flows_vec = common.molar_flows_vec
    z_vec = common.z_vec
    

    total_mass_flow_rate = 0.0u"kg/s"
    for (i, component) in enumerate(components)
        mw_i = clapeyron_model.params.Mw[i] * u"g/mol" |> u"kg/mol" 
        F_i = molar_flows_vec[i]
        total_mass_flow_rate += F_i * mw_i
    end

    avg_mw = total_mass_flow_rate / total_molar_flow

    superficial_mass_velocity = total_mass_flow_rate / reactor_cross_sectional_area
    
    molar_dens = molar_density(clapeyron_model, P, T, z_vec) 
    gas_density = molar_dens * avg_mw

    gas_viscosity = 2.0e-5u"Pa*s" 

    catalyst_particle_diameter = uconvert(u"m", catalyst_particle_diameter)

    term1 = 150 * (1 - bed_void_fraction)^2 / bed_void_fraction^3 * gas_viscosity * superficial_mass_velocity / (gas_density * catalyst_particle_diameter^2)
    term2 = 1.75 * (1 - bed_void_fraction) / bed_void_fraction^3 * superficial_mass_velocity^2 / (catalyst_particle_diameter * gas_density)
    dP_dL = -(term1 + term2) # This is in Pa/m
    #println("term1: ", uconvert(u"bar/m", term1))
    #println("term2: ", uconvert(u"bar/m", term2))
    #println("num: ", 1.75 * (1 - bed_void_fraction) * superficial_mass_velocity^2)
    #println("bed void fraction: ", bed_void_fraction)
    #println("superficial_mass_velocity: ", superficial_mass_velocity)
    #println("denom: ", (catalyst_particle_diameter * gas_density * bed_void_fraction^3))
    #println("catalyst_particle_diameter: ", catalyst_particle_diameter)
    #println("gas_density: ", gas_density)
    #println("bed_void_fraction: ", bed_void_fraction)
    println(println("dP_dL: ", uconvert(u"bar/m", dP_dL)))

    bed_density = (catalyst_density * (1 - bed_void_fraction))

    dP_dW = dP_dL / (bed_density * reactor_cross_sectional_area) # Pa/kg
    println("dP_dW: ", uconvert(u"bar/kg", dP_dW))
    return dP_dW
end



"""
I wonder for modelling isothermal, adiabatic, based on heat input, and pressure drop or no pressure drop should be modelled 
through multiple functions or one funciton with if statements
""" 
function pbr_ode_system!(du, u, p::PBRModel{Isothermal, ErgunPressureDrop}, catalyst_weight)
    T = p.thermal_model.reactor_temperature
    P = u[1]
    println("Next Iter")
    println("T, ", T)
    println("P, ", P)
    molar_flows_vec = u[2:end]
    
    common = calculate_rates_and_properties(T, P, molar_flows_vec, p)

    for (i, chemical) in enumerate(p.components)
        du[i + 1] = common.total_species_rates[chemical]
    end

    du[1] = calculate_dP_dW(u, p, T, P, common)
    return nothing
end

function pbr_ode_system!(du, u, p::PBRModel{NonIsothermalHX, ErgunPressureDrop}, catalyst_weight)
    T = u[1]
    P = u[2]
    println("Next Iter")
    println("T, ", T)
    println("P, ", P)
    molar_flows_vec = u[3:end]
    
    common = calculate_rates_and_properties(T, P, molar_flows_vec, p)

    for (i, chemical) in enumerate(p.components)
        du[i + 2] = common.total_species_rates[chemical]
    end

    du[1] = calculate_dT_dW(u, p, T, P, common)

    du[2] = calculate_dP_dW(u, p, T, P, common)
    return nothing
end

#The order in which reaction defs are put in the dictionary determines the order in which they will be applied

Wspan = (0.0u"kg", 100.0u"kg")

"""
heat_model = NonIsothermalHX(
    300u"K", #ambient temperature
    0.001u"W/(m^2*K)", #overall heat transfer coefficient
    0.001u"m^2/kg", #reactor surface area per kg of catalyst
    (6800u"W" / 100u"kg"), #heat wattage per kg
)
"""
heat_model = Isothermal(
    573.28u"K"
)

pressure_model = ErgunPressureDrop(
    5.0u"mm", #catalyst particle diameter
    0.3, #bed void fraction
    1000u"kg/m^3", #catalyst density 
    1.0u"m^2", #reactor cross sectional area
)

SMR_parameters = PBRModel(
    heat_model, #heat solver type
    pressure_model, #pressure solver type
    [MD_reaction, WGS_reaction], #reaction dict
    SMR_model, #thermo model
    ["CH3OH", "H2O", "CO", "H2", "CO2"], #components dict
    Dict("MD" => 90.7u"kJ/mol", "WGS" => -41.1u"kJ/mol"), #heat of reactions dict 
)

#u0 = vcat(reactor_inlet_temperature, reactor_inlet_pressure, SMR_molar_flows)
u0 = vcat(reactor_inlet_pressure, SMR_molar_flows)

pbr_ode_problem = ODEProblem(pbr_ode_system!, u0, Wspan, SMR_parameters)

pbr_sol = DifferentialEquations.solve(pbr_ode_problem, Tsit5())

temperature_over_reactor = [u[1] for u in pbr_sol.u]
pressure_over_reactor = [u[2] for u in pbr_sol.u]

mole_flow_Methanol_values = [u[3] for u in pbr_sol.u]
mole_flow_H2O_values = [u[4] for u in pbr_sol.u]
mole_flow_CO_values = [u[5] for u in pbr_sol.u]
mole_flow_H2_values = [u[6] for u in pbr_sol.u]
mole_flow_CO2_values = [u[7] for u in pbr_sol.u]

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
plot!(catalyst_weights, mole_flow_H2O_values, linewidth=2, marker=:circle, markersize=3, label="H2O")
plot!(catalyst_weights, mole_flow_CO_values, linewidth=2, marker=:circle, markersize=3, label="CO")
plot!(catalyst_weights, mole_flow_H2_values, linewidth=2, marker=:circle, markersize=3, label="H2")
plot!(catalyst_weights, mole_flow_CO2_values, linewidth=2, marker=:circle, markersize=3, label="CO2")


wattages_required = []
for i in eachindex(mole_flow_Methanol_values)
    push!(wattages_required, uconvert(u"W", (mole_flow_Methanol_values[1] - mole_flow_Methanol_values[i]) * 90.7u"kJ/mol"))
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

#display(plt)
    
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
plot!(catalyst_weights, mole_flow_H2O_values, linewidth=2, marker=:circle, markersize=3)
plot!(catalyst_weights, mole_flow_CO_values, linewidth=2, marker=:circle, markersize=3)
plot!(catalyst_weights, mole_flow_H2_values, linewidth=2, marker=:circle, markersize=3)
plot!(catalyst_weights, mole_flow_CO2_values, linewidth=2, marker=:circle, markersize=3)


plt2 = plot!(twinx(plt2), catalyst_weights, wattages_required,
    ylabel="Heat Input (W)",
    legend=:false, 
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    color=:red,
    framestyle=:box)

#display(plt2)

total_molar_flows_vec = []
concentrations_vec = []

pbr_sol.u[1]

for i in eachindex(pbr_sol.u)
        push!(total_molar_flows_vec, sum(pbr_sol.u[i][3:end]))
end

total_molar_flows_vec

for i in eachindex(pbr_sol.u)
    push!(concentrations_vec, [])
    for j in eachindex(pbr_sol.u[i][3:end])
        push!(concentrations_vec[i], pbr_sol.u[i][j+2] / total_molar_flows_vec[i])
    end
end

concentrations_vec

concentration_Methanol_values = [concentrations[1] for concentrations in concentrations_vec]
concentration_H2O_values = [concentrations[2] for concentrations in concentrations_vec]
concentration_CO_values = [concentrations[3] for concentrations in concentrations_vec]
concentration_H2_values = [concentrations[4] for concentrations in concentrations_vec]
concentration_CO2_values = [concentrations[5] for concentrations in concentrations_vec]


plt3 = plot(catalyst_weights, concentration_Methanol_values,
    xlabel="Catalyst Weight",
    ylabel="Molar Concentration Through Reactor",
    title="PBR Species Molar Concentration Through Reactor and \nHeat Input vs. Catalyst Weight",
    legend=:false, 
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    framestyle=:box)
    
# 2. Add the other species to the existing plot using plot!()
plot!(catalyst_weights, concentration_H2O_values, linewidth=2)
plot!(catalyst_weights, concentration_CO_values, linewidth=2)
plot!(catalyst_weights, concentration_H2_values, linewidth=2)
plot!(catalyst_weights, concentration_CO2_values, linewidth=2)

plt3 = plot!(twinx(plt3), catalyst_weights, wattages_required,
    ylabel="Heat Input (W)",
    legend=:false, 
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    color=:red,
    framestyle=:box)

#display(plt2)




#If you're getting bad conversion, it seems like the most probable cause is bad K_eq data for SMR

methanol_conversions_over_catalyst = []
for i in eachindex(pbr_sol.u)
    push!(methanol_conversions_over_catalyst, ((pbr_sol.u[1][3] - pbr_sol.u[i][3]) / pbr_sol.u[1][3]) * 100)
end

methanol_conversions_over_catalyst

plt = plot(catalyst_weights, methanol_conversions_over_catalyst,
    xlabel="Catalyst Weight",
    ylabel="Methanol Conversion %",
    title="PBR Methanol Conversion vs. Catalyst Weight",
    legend=:best, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    label="Methanol Conversion %",
    framestyle=:box)

  
water_conversions_over_catalyst = []
for i in eachindex(pbr_sol.u)
    push!(water_conversions_over_catalyst, ((pbr_sol.u[1][4] - pbr_sol.u[i][4]) / pbr_sol.u[1][4]) * 100)
end

water_conversions_over_catalyst
"""
plt = plot(catalyst_weights, water_conversions_over_catalyst,
    xlabel="Catalyst Weight",
    ylabel="Water Conversion %",
    title="PBR water Conversion vs. Catalyst Weight",
    legend=:best, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    label="Water Conversion %",
    framestyle=:box)
"""
