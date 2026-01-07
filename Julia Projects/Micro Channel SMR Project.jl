using Unitful, DifferentialEquations
using ComponentArrays
using Interpolations
using CoolProp
using Plots
using Roots

module GasProperties
    using Unitful

    # Standard Gas Constant
    const R_GAS = 8.314u"J/(mol*K)"

    # Shomate Equation Coefficients (NIST Webbook format) or NASA 7
    # Cp(t) = A + B*t + C*t^2 + D*t^3 + E/t^2 (J/mol*K) where t = T/1000
    # Valid roughly 298K - 1200K+
    const COEFFS = Dict(
        "methane" => (A=-0.703029, B=108.4773, C=-42.52157, D=5.862788, E=0.678515, MW=16.04u"g/mol"),
        "water"   => (A=30.09200, B=6.832514, C=6.793435, D=-2.534480, E=0.082139, MW=18.015u"g/mol"),
        "CO"      => (A=25.56759, B=6.096130, C=4.054656, D=-2.671301, E=0.131021, MW=28.01u"g/mol"),
        "hydrogen"=> (A=33.066178, B=-11.363417, C=11.432816, D=-2.772874, E=-0.158558, MW=2.016u"g/mol"),
        "CO2"     => (A=24.99735, B=55.18696, C=-33.69137, D=7.948387, E=-0.136638, MW=44.01u"g/mol"),
        "oxygen"  => (A=31.32234, B=-20.23531, C=57.86644, D=-36.50624, E=-0.007374, MW=31.99u"g/mol"),
        "nitrogen"=> (A=28.98641, B=1.853978, C=-9.647459, D=16.63537, E=0.000117, MW=28.014u"g/mol")
    )

    function get_cp_molar(species::String, T)
        t = ustrip(u"K", T) / 1000.0
        c = COEFFS[species]
        cp_val = c.A + c.B*t + c.C*t^2 + c.D*t^3 + c.E/t^2
        return cp_val * 1u"J/(mol*K)"
    end

    function get_mixture_cp(species_list, mole_fractions, T)
        cp_mix = 0.0u"J/(mol*K)"
        for (name, x) in zip(species_list, mole_fractions)
            if haskey(COEFFS, name)
                cp_mix += x * get_cp_molar(name, T)
            else
                # Fallback for solids like Carbon if ignored in fluid props
                cp_mix += x * 20.0u"J/(mol*K)" 
            end
        end
        return cp_mix
    end

    function get_mixture_density(species_list, mole_fractions, T, P)
        # Ideal Gas Law: rho_molar = P / (RT)
        # rho_mass = rho_molar * MW_mix
        
        avg_MW = 0.0u"g/mol"
        for (name, x) in zip(species_list, mole_fractions)
            if haskey(COEFFS, name)
                avg_MW += x * COEFFS[name].MW
            end
        end
        
        rho_molar = P / (R_GAS * T)
        return rho_molar * avg_MW
    end
    
    function get_mixture_molar_density(T, P)
        return P / (R_GAS * T)
    end
end

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
        forward_term *= pressure^nu
    end

    # 3. Calculate the reverse reaction "propensity" based on product partial pressures
    reverse_term = 1.0u"bar^0" # Initialize with a unitless 1
    for (product, nu) in reaction.products
        pressure = max(partial_pressures_dict[product], 0.0u"bar")
        reverse_term *= pressure^nu
    end
    
    # 4. Calculate the net rate of reaction (mol_rxn / kg_cat / s)
    net_reaction_rate = ((kf * forward_term) - (kr * reverse_term))
    
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
    end
    
    return [rates, net_reaction_rate]
end

function coolprop_mixture(species, mole_fractions)
    total_fraction = sum(mole_fractions)
    corrected_mole_fractions = [mole_fraction / total_fraction for mole_fraction in mole_fractions]

    coolprop_mixture_string = ""
    for i in eachindex(species)
        if i != length(species)
            coolprop_mixture_string *= species[i] * "[" * string(corrected_mole_fractions[i]) * "]" * "&"
        else 
            coolprop_mixture_string *= species[i] * "[" * string(corrected_mole_fractions[i]) * "]"
        end
    end
    return coolprop_mixture_string
end

function reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, characteristic_length)
    fluid_density = PropsSI("D", "T", fluid_temp, "P", fluid_pressure, fluid)
    dynamic_viscosity = PropsSI("V", "T", fluid_temp, "P", fluid_pressure, fluid)
    return (fluid_density * fluid_velocity * characteristic_length) / (dynamic_viscosity)
end

#overall_rate(SMR_reaction, SMR_model, reactor_inlet_temperature, 100000u"Pa", SMR_initial_C)

abstract type ThermalModel end
abstract type PressureModel end


struct Isothermal <: ThermalModel
    reactor_temperature::typeof(1.0u"K")
end

struct MixedThermalModel <: ThermalModel
    input_wattage_per_m::typeof(1.0u"W/m")

    ambient_temperature::typeof(1.0u"K")
    overall_heat_transfer_coeff::typeof(1.0u"W/(m^2*K)")
    reactor_surface_area_per_m::typeof(1.0u"m^2/m")
end

struct DarcyWeisbachPressureDrop <: PressureModel
    Re2f_interpolation::Interpolations.Extrapolation
    channel_hydraulic_diameter::typeof(1.0u"m")
    channel_cross_sectional_area::typeof(1.0u"m^2")
end

struct PBRModel{T<:ThermalModel, P<:PressureModel}
    thermal_model::T
    pressure_model::P
    reactions::Vector{HeterogeneousChemicalReaction}
    components::Vector{String}
    coolprop_components::Vector{String}
    heat_of_reactions_dict::Dict{String, typeof(1.0u"kJ/mol")}
end

function calculate_rates_and_properties(T::typeof(1.0u"K"), P, molar_flows_vec, p::PBRModel)
    reactions = p.reactions
    components = p.components

    total_molar_flow = sum(molar_flows_vec)
    molar_flows_dict = Dict(zip(components, molar_flows_vec))
    z_vec = [molar_flow / total_molar_flow for molar_flow in molar_flows_vec]
    
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

        for (chemical, rate) in species_rates
            total_species_rates[chemical] += rate
        end
    end
    return (net_reaction_rates=net_reaction_rates, total_species_rates=total_species_rates, molar_flows_dict=molar_flows_dict, molar_flows_vec=molar_flows_vec, z_vec=z_vec)
end

function calculate_dT_dL(p::PBRModel{MixedThermalModel}, T::typeof(1.0u"K"), P, common, dq_dL)
    reactions = p.reactions
    heat_of_reactions_dict = p.heat_of_reactions_dict

    molar_flows_dict = common.molar_flows_dict
    net_reaction_rates = common.net_reaction_rates

    #function specific variables
    input_wattage_per_m = p.thermal_model.input_wattage_per_m
    
    U = p.thermal_model.overall_heat_transfer_coeff
    a = p.thermal_model.reactor_surface_area_per_m
    ambient_temperature = p.thermal_model.ambient_temperature

    heat_exchange_with_environment = U*a*(ambient_temperature - T)

    catalyst_linear_density = 1.0u"kg/m" 

    delta_H_from_reactions = 0.0u"W/m"
    for i in eachindex(reactions)
        delta_H_from_reactions += heat_of_reactions_dict[reactions[i].name] * -net_reaction_rates[reactions[i].name] * catalyst_linear_density
    end

    coolprop_mixture_str = coolprop_mixture(p.coolprop_components, common.z_vec)

    fluid_molar_heat_capacity = 50u"J/(mol*K)"
    #PropsSI("CPMOLAR", "T", T, "P", P, coolprop_mixture_str)
    #This breaks with error message: CoolProp: One stationary point (not good) for 
    #T=573.28,p=2e+06,z=[ 0.847814688427, 0.127918202434, 0.0130353979437, 0.00331008440555, 0.00792162678894 ] : 
    #PropsSI("CPMOLAR","T",573.28,"P",2000000,"methane[0.24813895781637718]&water[0.7444168734491315]&CO[0.0024813895781637717]&hydrogen[0.0024813895781637717]&CO2[0.0024813895781637717]&Carbon_Solid[0.0]")
    
    total_heat_capacity_flow = sum([fluid_molar_heat_capacity * molar_flow for molar_flow in common.molar_flows_vec])
    
    dT_dL = (heat_exchange_with_environment + delta_H_from_reactions + input_wattage_per_m + dq_dL) / total_heat_capacity_flow
    
    return dT_dL 
end

function calculate_dP_dL(p::PBRModel{MixedThermalModel, DarcyWeisbachPressureDrop}, T::typeof(1.0u"K"), P, common)
    Re2f = p.pressure_model.Re2f_interpolation
    channel_hydraulic_diameter = p.pressure_model.channel_hydraulic_diameter
    channel_cross_sectional_area = p.pressure_model.channel_cross_sectional_area

    #pretty fragile but it works for now to get rid of carbon_solid 
    if length(common.z_vec) == 6
        coolprop_mixture_str = "SRK::" * coolprop_mixture(p.coolprop_components[1:2], common.z_vec[1:2])
    else
        coolprop_mixture_str = "SRK::" * coolprop_mixture(p.coolprop_components[1:2], common.z_vec[1:2])
    end

    #println(PropsSI("Tcrit", coolprop_mixture_str))
    #println(PropsSI("Pcrit", coolprop_mixture_str))
    #println(PropsSI("Dcrit", coolprop_mixture_str))
    fluid_molar_density = PropsSI("DMOLAR", "T", T, "P", P, coolprop_mixture_str)
    fluid_density = PropsSI("Dmass", "T", T, "P", P, coolprop_mixture_str)

    total_molar_flow = sum(common.molar_flows_vec)
    
    total_volumetric_flow = total_molar_flow / fluid_molar_density

    fluid_velocity = total_volumetric_flow / channel_cross_sectional_area

    #dynamic_viscosity = PropsSI("V", "T", T, "P", P, coolprop_mixture_str) #not avaliable for our mixture :(
    dynamic_viscosity = 3.5e-5u"Pa/s" #just using something close

    Re = ustrip((fluid_density * fluid_velocity * channel_hydraulic_diameter) / (dynamic_viscosity))

    if Re < 2200
        f = 64/Re # Laminar approximation
    else
        f = Re2f(Re)
    end

    dP_dL = -f * (fluid_density * fluid_velocity^2) / (2 * channel_hydraulic_diameter)

    return dP_dL
end
#=
function pbr_ode_system!(du, u, p::Vector, reactor_length)
    # Unpack parameters: p[1] is SMR, p[2] is Combustion
    smr_model = p[1]
    comb_model = p[2]
    
    # --- SMR System ---
    du.SMR_u.T = du.SMR_u.T .* 0.0
    du.SMR_u.P = du.SMR_u.P .* 0.0
    du.SMR_u.molar_flows_vec = du.SMR_u.molar_flows_vec .* 0.0


    SMR_T = u.SMR_u.T
    SMR_P = u.SMR_u.P
    SMR_molar_flows = u.SMR_u.molar_flows_vec
    
    common_SMR = calculate_rates_and_properties(SMR_T, SMR_P, SMR_molar_flows, smr_model)

    # Catalyst density assumption for rate conversion (mol/s/kg -> mol/s/m)
    catalyst_linear_density = 1.0u"kg/m"

    # Molar flow derivatives (reaction rates)
    for (i, chemical) in enumerate(smr_model.components)
        # du = rate * catalyst_density
        du.SMR_u.molar_flows_vec[i] = common_SMR.total_species_rates[chemical] * catalyst_linear_density
    end

    du.SMR_u.T = calculate_dT_dL(smr_model, SMR_T, SMR_P, common_SMR)
    du.SMR_u.P = calculate_dP_dL(smr_model, SMR_T, SMR_P, common_SMR)

    # --- Combustion System ---
    du.combustion_u.T = du.combustion_u.T .* 0.0
    du.combustion_u.P = du.combustion_u.P .* 0.0
    du.combustion_u.molar_flows_vec = du.combustion_u.molar_flows_vec .* 0.0


    comb_T = u.combustion_u.T
    comb_P = u.combustion_u.P
    comb_molar_flows = u.combustion_u.molar_flows_vec
    
    common_comb = calculate_rates_and_properties(comb_T, comb_P, comb_molar_flows, comb_model)

    for (i, chemical) in enumerate(comb_model.components)
        du.combustion_u.molar_flows_vec[i] = common_comb.total_species_rates[chemical] * catalyst_linear_density
    end

    du.combustion_u.T = calculate_dT_dL(comb_model, comb_T, comb_P, common_comb)
    du.combustion_u.P = calculate_dP_dL(comb_model, comb_T, comb_P, common_comb)

    return nothing
end
=#


function reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, characteristic_length)
    fluid_density = PropsSI("D", "T", fluid_temp, "P", fluid_pressure, fluid)
    dynamic_viscosity = 3.5e-5u"Pa/s"
    #PropsSI("V", "T", fluid_temp, "P", fluid_pressure, fluid)
    return ustrip((fluid_density * fluid_velocity * characteristic_length) / (dynamic_viscosity))
end
    
function prandtl_number(fluid, fluid_temperature, fluid_pressure)
    dynamic_viscosity = 3.5e-5u"Pa/s"
    #PropsSI("V", "T", fluid_temperature, "P", fluid_pressure, fluid)
    thermal_cond = 0.034u"W/(m*K)"
    #PropsSI("conductivity", "T", fluid_temperature, "P", fluid_pressure, fluid)
    specific_heat = PropsSI("C", "T", fluid_temperature, "P", fluid_pressure, fluid)
    return ustrip((dynamic_viscosity * specific_heat) / (thermal_cond))
end

function laminar_entry_Baehr_Stephan(Re, Pr, z, Di)
    eff_len = max(z, 1e-5) 
    Gz = Di/eff_len*Re*Pr
    return (3.657/tanh(2.264*Gz^(-1/3.)+ 1.7*Gz^(-2/3.0))
            + 0.0499*Gz*tanh(1.0/Gz))/tanh(2.432*Pr^(1/6.0)*Gz^(-1/6.0))
end

function turbulent_Gnielinski(Re, Pr, fd)
    return (fd/8.0)*(Re - 1e3)*Pr/(1.0 + 12.7*(fd/8.0)^0.5*(Pr^(2/3.0) - 1.0))
end


mutable struct Layer
    thickness::typeof(1.0u"m")
    thermal_conductivity::typeof(1.0u"W/(m*K)")
end

function get_linear_heat_transfer_coefficient(
        z,
        
        pipe_1_hydraulic_diameter, 
        pipe_1_fluid_thermal_cond,
        Re_1, Pr_1, fd_1,

        pipe_2_hydraulic_diameter, 
        pipe_2_fluid_thermal_cond,
        Re_2, Pr_2, fd_2,

        contact_width,
        
        layers
    )

    #pipe_1_fluid_thermal_cond = PropsSI("conductivity", "T", pipe_1_fluid_temp, "P", pipe_1_fluid_pressure, pipe_1_fluid)
    if Re_1 <= 2200 
        pipe_1_nusselt_number = laminar_entry_Baehr_Stephan(Re_1, Pr_1, z, pipe_1_hydraulic_diameter)
    else
        pipe_1_nusselt_number = turbulent_Gnielinski(Re_1, Pr_1, fd_1)
    end

    #pipe_2_fluid_thermal_cond = PropsSI("conductivity", "T", pipe_2_fluid_temp, "P", pipe_2_fluid_pressure, pipe_2_fluid)
    if Re_2 <= 2200 
        pipe_2_nusselt_number = laminar_entry_Baehr_Stephan(Re_2, Pr_2, z, pipe_2_hydraulic_diameter)
    else
        pipe_2_nusselt_number = turbulent_Gnielinski(Re_2, Pr_2, fd_2)
    end
    
    pipe_1_heat_transfer_coeff = (pipe_1_nusselt_number * pipe_1_fluid_thermal_cond) / pipe_1_hydraulic_diameter
    pipe_2_heat_transfer_coeff = (pipe_2_nusselt_number * pipe_2_fluid_thermal_cond) / pipe_2_hydraulic_diameter

    conv_resistance_1 = 1 / (pipe_1_heat_transfer_coeff * pipe_1_hydraulic_diameter)
    conv_resistance_2 = 1 / (pipe_2_heat_transfer_coeff * pipe_2_hydraulic_diameter)

    total_cond_resistance = sum([layer.thickness / (layer.thermal_conductivity * contact_width) for layer in layers])

    UA_prime = 1 / (conv_resistance_1 + total_cond_resistance + conv_resistance_2)
    
    return UA_prime
end

function calculate_UA_prime_system(
        layers, 
        model_1, model_1_common, model_1_T, model_1_P, model_1_molar_flows, 
        model_2, model_2_common, model_2_T, model_2_P, model_2_molar_flows, 
        z_position
    )
    # -- Pipe 1 -- 
    model_1_coolprop_mixture_str = "SRK::" * coolprop_mixture(model_1.coolprop_components, model_1_common.z_vec)
    fluid_1_thermal_conductivity = 0.034u"W/(m*K)"
    #PropsSI("conductivity", "T", model_1_T, "P", model_1_P, model_1_coolprop_mixture_str)

    fluid_1_molar_density = PropsSI("Dmolar", "T", model_1_T, "P", model_1_P, model_1_coolprop_mixture_str)
    fluid_1_density = PropsSI("Dmass", "T", model_1_T, "P", model_1_P, model_1_coolprop_mixture_str)

    fluid_1_total_molar_flow = sum(model_1_common.molar_flows_vec)
    
    fluid_1_total_volumetric_flow = fluid_1_total_molar_flow / fluid_1_molar_density

    fluid_velocity = fluid_1_total_volumetric_flow / model_1.pressure_model.channel_cross_sectional_area

    dynamic_viscosity = 3.5e-5u"Pa/s" #just using something close

    Re_1 = reynolds_number(model_1_coolprop_mixture_str, model_1_T, model_1_P, fluid_velocity, model_1.pressure_model.channel_hydraulic_diameter)
    Pr_1 = prandtl_number(model_1_coolprop_mixture_str, model_1_T, model_1_P)
    
    Re2f_1 = model_1.pressure_model.Re2f_interpolation
    if Re_1 < 2200
        fd_1 = 64/Re_1 # Laminar approximation
    else
        fd_1 = Re2f_1(Re_1)
    end

    # -- Pipe 2 -- 
    model_2_coolprop_mixture_str = "SRK::" * coolprop_mixture(model_2.coolprop_components, model_2_common.z_vec)
    fluid_2_thermal_conductivity = 0.034u"W/(m*K)"
    #PropsSI("conductivity", "T", model_2_T, "P", model_2_P, model_2_coolprop_mixture_str)

    fluid_2_molar_density = PropsSI("Dmolar", "T", model_2_T, "P", model_2_P, model_2_coolprop_mixture_str)
    fluid_2_density = PropsSI("Dmass", "T", model_2_T, "P", model_2_P, model_2_coolprop_mixture_str)

    fluid_2_total_molar_flow = sum(model_2_common.molar_flows_vec)
    
    fluid_2_total_volumetric_flow = fluid_2_total_molar_flow / fluid_2_molar_density

    fluid_velocity = fluid_2_total_volumetric_flow / model_2.pressure_model.channel_cross_sectional_area

    dynamic_viscosity = 3.5e-5u"Pa/s" #just using something close

    Re_2 = reynolds_number(model_2_coolprop_mixture_str, model_2_T, model_2_P, fluid_velocity, model_2.pressure_model.channel_hydraulic_diameter)
    Pr_2 = prandtl_number(model_2_coolprop_mixture_str, model_2_T, model_2_P)
    
    Re2f_2 = model_2.pressure_model.Re2f_interpolation
    if Re_2 < 2200
        fd_2 = 64/Re_2 # Laminar approximation
    else
        fd_2 = Re2f_2(Re_2)
    end

    phys_diameter = model_1.pressure_model.channel_hydraulic_diameter * (pi + 2) / pi

    return get_linear_heat_transfer_coefficient(
        z_position, 
        model_1.pressure_model.channel_hydraulic_diameter, fluid_1_thermal_conductivity, Re_1, Pr_1, fd_1, 
        model_2.pressure_model.channel_hydraulic_diameter, fluid_2_thermal_conductivity, Re_2, Pr_2, fd_2,
        phys_diameter,
        layers
    )
end

copper_layer = Layer(5u"mm", 401u"W/(m*K)")

layers = [copper_layer]

function pbr_ode_system!(du, u, p::Vector, z)
    # Unpack parameters: p[1] is SMR, p[2] is Combustion
    smr_model = p[1]
    comb_model = p[2]

    SMR_T = u[1]
    SMR_P = u[2]
    SMR_molar_flows = u[3:8]
    
    common_SMR = calculate_rates_and_properties(SMR_T, SMR_P, SMR_molar_flows, smr_model)

    # Catalyst density assumption for rate conversion (mol/s/kg -> mol/s/m)
    catalyst_linear_density = 1.0u"kg/m"

    for (i, chemical) in enumerate(smr_model.components)
        du[i+2] = common_SMR.total_species_rates[chemical] * catalyst_linear_density
    end
    
    du[2] = du[2] = calculate_dP_dL(smr_model, SMR_T, SMR_P, common_SMR)

    # --- Combustion System ---
    comb_T = u[9]
    comb_P = u[10]
    comb_molar_flows = u[11:end]
    
    common_comb = calculate_rates_and_properties(comb_T, comb_P, comb_molar_flows, comb_model)

    for (i, chemical) in enumerate(comb_model.components)
        du[i+10] = common_comb.total_species_rates[chemical] * catalyst_linear_density
    end
    
    du[10] = du[9] = calculate_dP_dL(comb_model, comb_T, comb_P, common_comb)

    #Calculate Heat Transfer due to difference in temperature between top and bottom reaction system
    layers = p[3]
    UA_prime = calculate_UA_prime_system(layers, 
                                    smr_model, common_SMR, SMR_T, SMR_P, SMR_molar_flows, 
                                    comb_model, common_comb, comb_T, comb_P, comb_molar_flows, 
                                    z
    )
    dq_dL = UA_prime * (comb_T - SMR_T)

    du[1] = calculate_dT_dL(smr_model, SMR_T, SMR_P, common_SMR, dq_dL)
    du[9] = calculate_dT_dL(comb_model, comb_T, comb_P, common_comb, dq_dL)
    return nothing
end

#The order in which reaction defs a put in the dictionary determines the order in which they will be applied

total_reactor_methane_mass_flow = 0.1u"g/s" 
SMR_to_combustor_methane_ratio = 3.0
SMR_methane_ratio = 1 / (SMR_to_combustor_methane_ratio + 1) * SMR_to_combustor_methane_ratio
combustor_methane_ratio = 1 / (SMR_to_combustor_methane_ratio + 1)

SMR_methane_molar_flow = (total_reactor_methane_mass_flow / 16.04u"g/mol") * SMR_methane_ratio
combustor_methane_molar_flow = (total_reactor_methane_mass_flow / 16.04u"g/mol") * combustor_methane_ratio

# -- START OF SMR Definition --

#chemistry parameters
reactor_inlet_temperature = 300.13u"°C" |> u"K"
reactor_inlet_pressure = 1.0u"bar" # SMR usually runs at higher pressure (15-30 bar), updated from 1 bar

# Base molar flow reference

#0.0652u"mol/s"

# Inlet Composition: Steam to Carbon ratio ~ 3:1 is typical to prevent coking
# Order: [CH4, H2O, CO, H2, CO2, Cs]
# Cs (Solid Carbon) is 0 at inlet
z_smr = [1.0, 3.0, 0.01, 0.01, 0.01, 0.0] 
SMR_initial_C = z_smr .* u"mol/L"

# --- Reaction 1: Overall SMR ---
# CH4 + H2O <-> CO + 3H2
ref_T_SMR = 700.0u"°C" |> u"K" # Ref T higher for SMR parameters
# Est. Kinetic parameters
kf_ref_smr = 5.0e-1u"mol/(kg*s*bar^2)" 
Ea_f_smr = 240.0u"kJ/mol" 
kf_A_smr = arrenhius_equation_pre_exponential_factor(kf_ref_smr, Ea_f_smr, ref_T_SMR)

# Thermodynamics for SMR
# dH ~ +206 kJ/mol (Endothermic)
# dG(298) ~ +142 kJ/mol
Keq_ref_smr = K_gibbs_free(298u"K", ref_T_SMR, 142.0u"kJ/mol", 206.0u"kJ/mol")

kr_ref_smr = kf_ref_smr / Keq_ref_smr
Ea_r_smr = Ea_f_smr - 206.0u"kJ/mol" # Ea_r approx
kr_A_smr = arrenhius_equation_pre_exponential_factor(kr_ref_smr, Ea_r_smr, ref_T_SMR) / 1.0u"bar^2"

SMR_reaction = HeterogeneousChemicalReaction(
    "SMR",
    Dict("CH4" => 1, "H2O" => 1),
    Dict("CO" => 1, "H2" => 3),
    ["CH4", "H2O", "CO", "H2"],
    kf_A_smr, Ea_f_smr,
    kr_A_smr, Ea_r_smr
)

# --- Reaction 2: Water Gas Shift (WGS) ---
# CO + H2O <-> CO2 + H2
ref_T_WGS = 300.13u"°C" |> u"K"
kf_ref_wgs = 0.2e-1u"mol/(kg*s*bar^2)"
Ea_f_wgs = 60u"kJ/mol"
kf_A_wgs = arrenhius_equation_pre_exponential_factor(kf_ref_wgs, Ea_f_wgs, ref_T_WGS)

Keq_ref_wgs = K_gibbs_free(298u"K", ref_T_WGS, -28.63u"kJ/mol", -41.1u"kJ/mol")

kr_ref_wgs = kf_ref_wgs / Keq_ref_wgs
Ea_r_wgs = 100.0u"kJ/mol" # or Ea_f - dH
kr_A_wgs = arrenhius_equation_pre_exponential_factor(kr_ref_wgs, Ea_r_wgs, ref_T_WGS)

WGS_reaction = HeterogeneousChemicalReaction(
    "WGS",
    Dict("CO" => 1, "H2O" => 1),
    Dict("CO2" => 1, "H2" => 1),
    ["CO", "H2O", "CO2", "H2"],
    kf_A_wgs, Ea_f_wgs,
    kr_A_wgs, Ea_r_wgs
)

# --- Reaction 3: Coking (Methane Cracking) ---
# CH4 -> C(s) + 2H2
# Note: Irreversible approximation for this estimation
ref_T_Coking = 700.0u"°C" |> u"K"
kf_ref_coke = 1.0e-4u"mol/(kg*s*bar)" # Usually slow compared to SMR if steam is high
Ea_f_coke = 150.0u"kJ/mol"
kf_A_coke = arrenhius_equation_pre_exponential_factor(kf_ref_coke, Ea_f_coke, ref_T_Coking)

# Reverse parameters set to effectively zero for irreversible assumption
kr_A_coke = 0.0u"mol/(kg*s*bar^3)" 
Ea_r_coke = 0.0u"kJ/mol"

Coking_reaction = HeterogeneousChemicalReaction(
    "Coking",
    Dict("CH4" => 1),
    Dict("Cs" => 1, "H2" => 2),
    ["CH4", "Cs", "H2"],
    kf_A_coke, Ea_f_coke,
    kr_A_coke, Ea_r_coke
)

# Flow calculations
molar_flow_correction_factor_smr = SMR_methane_molar_flow / SMR_initial_C[1]
SMR_molar_flows = SMR_initial_C .* molar_flow_correction_factor_smr

# -- Physical Model Parameters --
input_wattage_per_m = 0u"W/m"
ambient_temperature = 25.13u"°C" |> u"K"
overall_heat_trans_coeff = 60u"W/(m^2*K)"
reactor_surface_area_per_m = 6u"cm^2/cm"

SMR_thermal_model = MixedThermalModel(
    input_wattage_per_m, 
    ambient_temperature,
    overall_heat_trans_coeff,
    reactor_surface_area_per_m
)

# Hydraulic / Pressure Drop Setup
Re_range = 2200:100:50000000
f_vals = []
for Re in Re_range
    find_f = x -> 2 * log10(Re*sqrt(x)) - 0.8 - 1/sqrt(x)
    f = find_zero(find_f, 1e-3)
    push!(f_vals, f)
end
Re2f = linear_interpolation(Re_range, f_vals)

channel_diameter = 3.0u"mm"
channel_radius = channel_diameter / 2
area = (pi * channel_radius^2) / 2
wetted_perimeter = 1.22 * channel_radius
channel_hydraulic_diameter = (4 * area) / wetted_perimeter
channel_cross_sectional_area = area / 2

SMR_pressure_model = DarcyWeisbachPressureDrop(
    Re2f,
    channel_hydraulic_diameter, 
    channel_cross_sectional_area
)

SMR_parameters = PBRModel(
    SMR_thermal_model, 
    SMR_pressure_model, 
    [SMR_reaction, WGS_reaction, Coking_reaction], 
    ["CH4", "H2O", "CO", "H2", "CO2", "Cs"], # Component Order
    ["methane", "water", "CO", "hydrogen", "CO2"], # Names
    Dict("SMR" => 206.0u"kJ/mol", "WGS" => -41.1u"kJ/mol", "Coking" => 74.9u"kJ/mol"), 
)

# -- END OF SMR Definition --


# -- START OF Combustion Definition --

combustion_inlet_temperature = 300.0u"°C" |> u"K"
combustion_inlet_pressure = 1.0u"bar"

# Combustion Inlet: Methane + Air (approx 21% O2, 79% N2)
# Stoichiometric: CH4 + 2O2 -> CO2 + 2H2O
# Using slight excess air (Lean) to ensure full combustion
z_comb = [1.0, 2.2, 8.3, 0.001, 0.001] # [CH4, O2, N2, CO2, H2O]
comb_initial_C = z_comb .* u"mol/L"
comb_correction = combustor_methane_molar_flow / comb_initial_C[1]
combustion_molar_flows = comb_initial_C .* comb_correction

# --- Reaction: Methane Combustion ---
# CH4 + 2O2 -> CO2 + 2H2O
ref_T_Comb = 500.0u"°C" |> u"K"
kf_ref_comb = 1.0e-1u"mol/(kg*s*bar^3)" 
Ea_f_comb = 80.0u"kJ/mol"
kf_A_comb = arrenhius_equation_pre_exponential_factor(kf_ref_comb, Ea_f_comb, ref_T_Comb)

# Combustion is effectively irreversible
kr_A_comb = 0.0u"mol/(kg*s*bar^3)"
Ea_r_comb = 0.0u"kJ/mol"

Combustion_Reaction = HeterogeneousChemicalReaction(
    "Combustion",
    Dict("CH4" => 1, "O2" => 2),
    Dict("CO2" => 1, "H2O" => 2),
    ["CH4", "O2", "CO2", "H2O"],
    kf_A_comb, Ea_f_comb,
    kr_A_comb, Ea_r_comb
)

combustion_thermal_model = MixedThermalModel(
    input_wattage_per_m, 
    ambient_temperature,
    overall_heat_trans_coeff,
    reactor_surface_area_per_m
)

combustion_pressure_model = DarcyWeisbachPressureDrop(
    Re2f,
    channel_hydraulic_diameter, 
    channel_cross_sectional_area
)

combustion_parameters = PBRModel(
    combustion_thermal_model, 
    combustion_pressure_model, 
    [Combustion_Reaction], 
    ["CH4", "O2", "N2", "CO2", "H2O"], # N2 is inert
    ["methane", "oxygen", "nitrogen", "CO2", "water"],
    Dict("Combustion" => -802.0u"kJ/mol"), 
)

# -- END OF Combustion Definition --

parameter_vec = [SMR_parameters, combustion_parameters, layers]

#SMR_u0 = ComponentVector(T = reactor_inlet_temperature, P = reactor_inlet_pressure, molar_flows_vec = SMR_molar_flows)

#combustion_u0 = ComponentVector(T = combustion_inlet_temperature, P = combustion_inlet_pressure, molar_flows_vec = combustion_molar_flows)

#PBR_u0 = ComponentVector(SMR_u = SMR_u0, combustion_u = combustion_u0)
u0 = vcat(reactor_inlet_temperature, reactor_inlet_pressure, SMR_molar_flows, combustion_inlet_temperature, combustion_inlet_pressure, combustion_molar_flows)

lspan = (0.0u"m", 1.0u"m")
pbr_ode_problem = ODEProblem(pbr_ode_system!, u0, lspan, parameter_vec)

pbr_sol = DifferentialEquations.solve(pbr_ode_problem, Tsit5())

#TODO
#sub reversible_pressure_rate_law for Langmuir-Hinshelwood-Hougen-Watson (LHHW) kinetics (e.g., the Xu-Froment model).