using Revise

using Unitful, DifferentialEquations
using Clapeyron
using Plots

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

#An isothermal dT_dW is not required because there's not temeprature change

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
    sub_heat_capacities = [90u"J/K", 70u"J/K", 20u"J/K", 20u"J/K", 20u"J/K"]
    for component in components
        current_model = PR([clapeyron_model.components[i]])
        try 
            comp_heat_capacity = isobaric_heat_capacity(current_model, P, T, [1.0]) #Not sure if I should use isobaric or isochoric
            println(comp_heat_capacity)
        catch  
            comp_heat_capacity = sub_heat_capacities[i]
        end
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

function calculate_dT_dW(u, p::PBRModel{InputWattage}, T::typeof(1.0u"K"), P, common)
    clapeyron_model = p.clapeyron_model
    components = p.components
    reactions = p.reactions
    heat_of_reactions_dict = p.heat_of_reactions_dict
    
    input_wattage_per_kg = p.thermal_model.input_wattage_per_kg

    molar_flows_dict = common.molar_flows_dict
    net_reaction_rates = common.net_reaction_rates
    
    delta_H_from_reactions = 0.0u"W/kg"
    for i in eachindex(reactions)
        println("current_heat_of_reaction: ", heat_of_reactions_dict[reactions[i].name])
        println("-net_reaction_rates: ", -net_reaction_rates[reactions[i].name])
        delta_H_from_reactions += heat_of_reactions_dict[reactions[i].name] * -net_reaction_rates[reactions[i].name]
    end

    total_heat_capacity_flow_rate = 0.0u"J/(K*s)"
    i = 1
    sub_heat_capacities = [90u"J/K", 70u"J/K", 20u"J/K", 20u"J/K", 20u"J/K"]
    for component in components
        current_model = PR([clapeyron_model.components[i]])
        try 
            comp_heat_capacity = isobaric_heat_capacity(current_model, P, T, [1.0]) #Not sure if I should use isobaric or isochoric
            println(comp_heat_capacity)
        catch  
            comp_heat_capacity = sub_heat_capacities[i]
        end
        comp_heat_capacity /= 1.0u"mol" #for some reason clapeyron doesn't return this in moles just j/k
        total_heat_capacity_flow_rate += molar_flows_dict[component] * comp_heat_capacity 
        i += 1
    end

    dT_dW = (delta_H_from_reactions + input_wattage_per_kg) / total_heat_capacity_flow_rate

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
    println("HELLO", typeof(u))

    total_mass_flow_rate = 0.0u"kg/s"
    for (i, component) in enumerate(components)
        mw_i = clapeyron_model.params.Mw[i] * u"g/mol" |> u"kg/mol" 
        F_i = molar_flows_vec[i]
        total_mass_flow_rate += F_i * mw_i
    end

    avg_mw = total_mass_flow_rate / total_molar_flow

    superficial_mass_velocity = total_mass_flow_rate / reactor_cross_sectional_area
    #try 
        #molar_dens = molar_density(clapeyron_model, P, T, z_vec)
    #finally
        
    #end
    molar_dens = 47.5u"mol/m^3"
    
    gas_density = molar_dens * avg_mw

    gas_viscosity = 2.0e-5u"Pa*s" 

    catalyst_particle_diameter = uconvert(u"m", catalyst_particle_diameter)

    term1 = 150 * (1 - bed_void_fraction)^2 / bed_void_fraction^3 * gas_viscosity * superficial_mass_velocity / (gas_density * catalyst_particle_diameter^2)
    term2 = 1.75 * (1 - bed_void_fraction) / bed_void_fraction^3 * superficial_mass_velocity^2 / (catalyst_particle_diameter * gas_density)
    dP_dL = -(term1 + term2) # This is in Pa/m

    bed_density = (catalyst_density * (1 - bed_void_fraction))

    dP_dW = dP_dL / (bed_density * reactor_cross_sectional_area) # Pa/kg
    return dP_dW
end

function calculate_dP_dW(u, p::PBRModel{NonIsothermal, ErgunPressureDrop}, T::typeof(1.0u"K"), P, common)
    clapeyron_model = p.clapeyron_model
    components = p.components
    reactor_cross_sectional_area = p.pressure_model.reactor_cross_sectional_area
    bed_void_fraction = p.pressure_model.bed_void_fraction
    catalyst_particle_diameter = p.pressure_model.catalyst_particle_diameter
    catalyst_density = p.pressure_model.catalyst_density 

    molar_flows_vec = common.molar_flows_vec
    z_vec = common.z_vec
    println("HELLO", typeof(u))

    total_mass_flow_rate = 0.0u"kg/s"
    for (i, component) in enumerate(components)
        mw_i = clapeyron_model.params.Mw[i] * u"g/mol" |> u"kg/mol" 
        F_i = molar_flows_vec[i]
        total_mass_flow_rate += F_i * mw_i
    end

    avg_mw = total_mass_flow_rate / total_molar_flow

    superficial_mass_velocity = total_mass_flow_rate / reactor_cross_sectional_area
    #try 
        #molar_dens = molar_density(clapeyron_model, P, T, z_vec)
    #finally
        
    #end
    molar_dens = 47.5u"mol/m^3"
    
    gas_density = molar_dens * avg_mw

    gas_viscosity = 2.0e-5u"Pa*s" 

    catalyst_particle_diameter = uconvert(u"m", catalyst_particle_diameter)

    term1 = 150 * (1 - bed_void_fraction)^2 / bed_void_fraction^3 * gas_viscosity * superficial_mass_velocity / (gas_density * catalyst_particle_diameter^2)
    term2 = 1.75 * (1 - bed_void_fraction) / bed_void_fraction^3 * superficial_mass_velocity^2 / (catalyst_particle_diameter * gas_density)
    dP_dL = -(term1 + term2) # This is in Pa/m

    bed_density = (catalyst_density * (1 - bed_void_fraction))

    dP_dW = dP_dL / (bed_density * reactor_cross_sectional_area) # Pa/kg
    return dP_dW
end

function calculate_dP_dW(u, p::PBRModel{InputWattage, ErgunPressureDrop}, T::typeof(1.0u"K"), P, common)
    clapeyron_model = p.clapeyron_model
    components = p.components
    reactor_cross_sectional_area = p.pressure_model.reactor_cross_sectional_area
    bed_void_fraction = p.pressure_model.bed_void_fraction
    catalyst_particle_diameter = p.pressure_model.catalyst_particle_diameter
    catalyst_density = p.pressure_model.catalyst_density 

    molar_flows_vec = common.molar_flows_vec
    z_vec = common.z_vec
    println("HELLO", typeof(u))

    total_mass_flow_rate = 0.0u"kg/s"
    for (i, component) in enumerate(components)
        mw_i = clapeyron_model.params.Mw[i] * u"g/mol" |> u"kg/mol" 
        F_i = molar_flows_vec[i]
        total_mass_flow_rate += F_i * mw_i
    end

    avg_mw = total_mass_flow_rate / total_molar_flow

    superficial_mass_velocity = total_mass_flow_rate / reactor_cross_sectional_area
    #try 
        #molar_dens = molar_density(clapeyron_model, P, T, z_vec)
    #finally
        
    #end
    molar_dens = 47.5u"mol/m^3"
    
    gas_density = molar_dens * avg_mw

    gas_viscosity = 2.0e-5u"Pa*s" 

    catalyst_particle_diameter = uconvert(u"m", catalyst_particle_diameter)

    term1 = 150 * (1 - bed_void_fraction)^2 / bed_void_fraction^3 * gas_viscosity * superficial_mass_velocity / (gas_density * catalyst_particle_diameter^2)
    term2 = 1.75 * (1 - bed_void_fraction) / bed_void_fraction^3 * superficial_mass_velocity^2 / (catalyst_particle_diameter * gas_density)
    dP_dL = -(term1 + term2) # This is in Pa/m

    bed_density = (catalyst_density * (1 - bed_void_fraction))

    dP_dW = dP_dL / (bed_density * reactor_cross_sectional_area) # Pa/kg
    return dP_dW
end


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

function pbr_ode_system!(du, u, p::PBRModel{NonIsothermal, ErgunPressureDrop}, catalyst_weight)
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

function pbr_ode_system!(du, u, p::PBRModel{InputWattage, ErgunPressureDrop}, catalyst_weight)
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

#The order in which reaction defs a put in the dictionary determines the order in which they will be applied

reactor_inlet_temperature = 350.13u"°C" |> u"K"

reactor_inlet_pressure = 1.0u"bar"

SMR_methanol_molar_flow = 0.06520168451976938u"mol/s"

SMR_model = PR(["Methanol", "water", "carbon monoxide", "hydrogen", "Carbon Dioxide"])

z = [1.0, 1.3, 0.001, 0.001, 0.001]
#z = [1, 1.3, 0.1, 0.1, 0.001]
SMR_initial_C = z .* u"mol/L"

# CH3OH -> CO + 2H2
# Below is for MD
ref_T = 300.13u"°C" |> u"K"
kf_ref = 1e-1u"mol/(kg*s*bar)" 
Ea_f = 90u"kJ/mol" 
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

Keq_ref = K_gibbs_free(298u"K", ref_T, 25.2u"kJ/mol", 90.7u"kJ/mol")
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
ref_T = 300.13u"°C" |> u"K"
kf_ref = 0.2e-1u"mol/(kg*s*bar^2)" #around 1e-6 to 1e-5
Ea_f = 60u"kJ/mol" #usually around 55-100 kJ/mol
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

# For reverse reaction, calculate kr from Keq. 
Keq_ref = K_gibbs_free(298u"K", ref_T, -28.63u"kJ/mol", -41.1u"kJ/mol")
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

molar_flow_correction_factor = SMR_methanol_molar_flow / SMR_initial_C[1]

SMR_molar_flows = SMR_initial_C .* molar_flow_correction_factor

total_molar_flow = sum(SMR_molar_flows)


Wspan = (0.0u"kg", 50u"kg")

thermal_model = Isothermal(
    reactor_inlet_temperature
)

pressure_model = ErgunPressureDrop(
    5.0u"mm", #catalyst particle diameter
    0.3, #bed void fraction
    1000u"kg/m^3", #catalyst density 
    1.0u"m^2", #reactor cross sectional area
)

SMR_parameters = PBRModel(
    thermal_model, #heat solver type
    pressure_model, #pressure solver type
    [MD_reaction, WGS_reaction], #reaction dict
    SMR_model, #thermo model
    ["CH3OH", "H2O", "CO", "H2", "CO2"], #components dict
    Dict("MD" => 90.7u"kJ/mol", "WGS" => -41.1u"kJ/mol"), #heat of reactions dict 
)

if typeof(SMR_parameters.thermal_model) == Isothermal
    u0 = vcat(reactor_inlet_pressure, SMR_molar_flows)
else
    u0 = vcat(reactor_inlet_temperature, reactor_inlet_pressure, SMR_molar_flows)
end

pbr_ode_problem = ODEProblem(pbr_ode_system!, u0, Wspan, SMR_parameters)

pbr_sol = DifferentialEquations.solve(pbr_ode_problem, Tsit5())

u_start = length(u0) - length(SMR_parameters.components) + 1

if typeof(SMR_parameters.thermal_model) == Isothermal
    pressure_over_reactor = [u[1] for u in pbr_sol.u]

    mole_flow_Methanol_values = [u[2] for u in pbr_sol.u]
    mole_flow_H2O_values = [u[3] for u in pbr_sol.u]
    mole_flow_CO_values = [u[4] for u in pbr_sol.u]
    mole_flow_H2_values = [u[5] for u in pbr_sol.u]
    mole_flow_CO2_values = [u[6] for u in pbr_sol.u]
else
    temperature_over_reactor = [u[1] for u in pbr_sol.u]
    pressure_over_reactor = [u[2] for u in pbr_sol.u]

    mole_flow_Methanol_values = [u[3] for u in pbr_sol.u]
    mole_flow_H2O_values = [u[4] for u in pbr_sol.u]
    mole_flow_CO_values = [u[5] for u in pbr_sol.u]
    mole_flow_H2_values = [u[6] for u in pbr_sol.u]
    mole_flow_CO2_values = [u[7] for u in pbr_sol.u]
end


mole_flow_Methanol_values
mole_flow_H2O_values
mole_flow_CO_values
mole_flow_H2_values
mole_flow_CO2_values

catalyst_weights = pbr_sol.t

packing_density = 0.5u"cm^3/cm^3"

catalyst_density = 8.96u"g/cm^3"

approximate_reactor_volues = catalyst_weights ./ (catalyst_density * packing_density)
approximate_reactor_volues = uconvert.(u"L", approximate_reactor_volues)

plt = plot(catalyst_weights, mole_flow_Methanol_values,
    xlabel="Catalyst Weight",
    ylabel="Molar Flow Through Reactor",
    title="PBR Species Molar Flow Through \n Reactor vs. Catalyst Weight at 350°C",
    legend=:best, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    label="CH3OH",
    framestyle=:box,
    dpi=1000)
    
# 2. Add the other species to the existing plot using plot!()
plot!(catalyst_weights, mole_flow_H2O_values, linewidth=2, marker=:circle, markersize=3, label="H2O")
plot!(catalyst_weights, mole_flow_CO_values, linewidth=2, marker=:circle, markersize=3, label="CO")
plot!(catalyst_weights, mole_flow_H2_values, linewidth=2, marker=:circle, markersize=3, label="H2")
plot!(catalyst_weights, mole_flow_CO2_values, linewidth=2, marker=:circle, markersize=3, label="CO2")

savefig(plt, "Math IA PBR 350C Long.png")

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

for i in eachindex(pbr_sol.u)
    push!(total_molar_flows_vec, sum(pbr_sol.u[i][u_start:end]))
end

total_molar_flows_vec

"""
mass_flow_consistentcy_check = [] #0.0u"kg/s"
for i in eachindex(pbr_sol.u)
    push!(mass_flow_consistentcy_check, sum(mole_flow_Methanol_values * 32.04u"g/mol" + mole_flow_H2O_values * 18.01528u"g/mol" + mole_flow_CO_values * 28.01u"g/mol" + mole_flow_H2_values * 2.016u"g/mol" + mole_flow_CO2_values * 44.009u"g/mol"))
end
mass_flow_consistentcy_check
"""

for i in eachindex(pbr_sol.u)
    push!(concentrations_vec, [])
    for j in eachindex(pbr_sol.u[i][u_start:end])
        push!(concentrations_vec[i], pbr_sol.u[i][j+u_start-1] / total_molar_flows_vec[i])
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
    push!(methanol_conversions_over_catalyst, ((mole_flow_Methanol_values[1] - mole_flow_Methanol_values[i]) / mole_flow_Methanol_values[1]) * 100)
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
    push!(water_conversions_over_catalyst, ((mole_flow_H2O_values[1] - mole_flow_H2O_values[i]) / mole_flow_H2O_values[1]) * 100)
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
#Just to show that CO2 Values are actually changing


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
    temperature = reactor_inlet_temperature,
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
        temperature = reactor_inlet_temperature, # Isothermal assumption for now
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
ref_T = reactor_inlet_temperature
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

#rate_law_for_chemicals(WGS_reaction, WGS_model, reactor_inlet_temperature, 100000u"Pa", WGS_initial_C)


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




"""
using NLsolve
using Symbolics
function isothermal_pbr_ode_system!(du, u, p, catalyst_weight)
    reactions = p.reactions
    model = p.clapeyron_model
    components = p.components
    heat_of_reactions_dict = p.heat_of_reactions_dict
    ambient_temperature = p.ambient_temperature
    U = p.overall_heat_transfer_coeff
    a = p.reactor_surface_area 
    reactor_cross_sectional_area = p.reactor_cross_sectional_area
    catalyst_particle_diameter = p.catalyst_particle_diameter
    bed_void_fraction = p.bed_void_fraction
    catalyst_density = p.catalyst_density

    T = p.reactor_temperature
    P = u[1]
    println("Next Iter")
    println("T, ", T)
    println("P, ", P)
    molar_flows_vec = u[2:end]

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

    for (i, chemical) in enumerate(components)
        du[i + 1] = total_species_rates[chemical]
    end

    total_mass_flow_rate = 0.0u"kg/s"
    for (i, component) in enumerate(components)
        mw_i = model.params.Mw[i] * u"g/mol" |> u"kg/mol" 
        F_i = molar_flows_vec[i]
        total_mass_flow_rate += F_i * mw_i
    end

    avg_mw = total_mass_flow_rate / total_molar_flow

    superficial_mass_velocity = total_mass_flow_rate / reactor_cross_sectional_area
    
    molar_dens = molar_density(model, P, T, z_vec) 
    gas_density = molar_dens * avg_mw

    gas_viscosity = 2.0e-5u"Pa*s" 

    catalyst_particle_diameter = uconvert(u"m", catalyst_particle_diameter)

    term1 = 150 * (1 - bed_void_fraction)^2 / bed_void_fraction^3 * gas_viscosity * superficial_mass_velocity / (gas_density * catalyst_particle_diameter^2)
    term2 = 1.75 * (1 - bed_void_fraction) / bed_void_fraction^3 * superficial_mass_velocity^2 / (catalyst_particle_diameter * gas_density)
    dP_dL = -(term1 + term2) # This is in Pa/m

    bed_density = (catalyst_density * (1 - bed_void_fraction))

    dP_dW = dP_dL / (bed_density * reactor_cross_sectional_area) # Pa/kg
    println("dP_dW: ", uconvert(u"bar/kg", dP_dW))
    println()
    println()

    du[1] = dP_dW
    return nothing
end

reactor_inlet_temperature = 1000u"K"

#The order in which reaction defs are put in the dictionary determines the order in which they will be applied
SMR_parameters = (
    reactions = [SMR_reaction, WGS_reaction], 
    clapeyron_model = SMR_model, 
    components = ["CH3OH", "H2O", "CO", "H2", "CO2"],
    heat_of_reactions_dict = Dict("SMR" => 49.4u"kJ/mol", "WGS" => -41.1u"kJ/mol"),
    ambient_temperature = 300u"K",
    overall_heat_transfer_coeff = 0.001u"W/(m^2*K)",
    reactor_surface_area = 0.001u"m^2/kg",
    reactor_cross_sectional_area = 1.0u"m^2",
    catalyst_particle_diameter = 5.0u"mm",
    bed_void_fraction = 0.3,
    catalyst_density = 1000u"kg/m^3",
    reactor_temperature = reactor_inlet_temperature,
)

Wspan = (0.0u"kg", 0.01u"kg")

u0 = vcat(reactor_inlet_pressure, SMR_molar_flows)

pbr_ode_problem = ODEProblem(isothermal_pbr_ode_system!, u0, Wspan, SMR_parameters)

pbr_sol = DifferentialEquations.solve(pbr_ode_problem, Tsit5())

pressure_over_reactor = [u[1] for u in pbr_sol.u]

mole_flow_Methanol_values = [u[2] for u in pbr_sol.u]
mole_flow_H2O_values = [u[3] for u in pbr_sol.u]
mole_flow_CO_values = [u[4] for u in pbr_sol.u]
mole_flow_H2_values = [u[5] for u in pbr_sol.u]
mole_flow_CO2_values = [u[6] for u in pbr_sol.u]

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

#If you're getting bad conversion, it seems like the most probable cause is bad K_eq data for SMR

methanol_conversions_over_catalyst = []
for i in eachindex(pbr_sol.u)
    push!(methanol_conversions_over_catalyst, ((pbr_sol.u[1][2] - pbr_sol.u[i][2]) / pbr_sol.u[1][3]) * 100)
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
#Just to show that CO2 Values are actually changing
"""