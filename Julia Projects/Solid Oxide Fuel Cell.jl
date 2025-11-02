using Unitful
#using Symbolics, Nemo, Groebner
using Clapeyron
#using NLsolve

#using CoolProp : cp
import CoolProp as Cp

desired_final_DC_wattage = 25000u"W"

desired_final_wattage = desired_final_DC_wattage / 0.97 #Account for AC losses

stack_thermal_efficiency = 0.55

stack_total_wattage = desired_final_wattage / stack_thermal_efficiency

stack_chemical_waste_heat = stack_total_wattage - desired_final_wattage 

desired_final_voltage = 240u"V"

desired_cell_amperage = uconvert(u"A", desired_final_wattage / desired_final_voltage)

function Nerst_eq(standard_potential_at_operating_temperature, stack_temperature,
                    number_of_moles_of_electrons_transferred_per_mole_of_reaction,
                    water_molar_concentration, hydrogen_molar_concentration, oxygen_molar_concentration)
    R = 8.314u"J/(mol*K)"
    standard_potential_at_operating_temperature - ((R * stack_temperature) / (number_of_moles_of_electrons_transferred_per_mole_of_reaction * faraday_constant)) * log((water_molar_concentration / (hydrogen_molar_concentration * oxygen_molar_concentration^0.5)))
end

faraday_constant = 96485u"C/mol"
theoretical_cell_voltage = uconvert(u"V", Nerst_eq(1.03u"V", (750u"°C" |> u"K"), 2, 0.15, 0.5, 0.21))

activation_voltage_loss = 0.15u"V" #0.1-0.25
ohmic_voltage_loss = 0.10u"V" #0.05-0.15
concentration_voltage_loss = 0.04u"V" #0.02-0.08

total_losses = activation_voltage_loss + ohmic_voltage_loss + concentration_voltage_loss

actual_cell_voltage = theoretical_cell_voltage - total_losses #about 1 volt

required_number_of_cells = ceil(desired_final_voltage / actual_cell_voltage) #Apply a ceiling because you can't have a fraction of a cell

final_cell_amperage = uconvert(u"A", (desired_final_wattage / actual_cell_voltage) / required_number_of_cells)

resistive_heat = uconvert(u"W", (final_cell_amperage^2) * ((total_losses * required_number_of_cells) / final_cell_amperage)) 
#Note that the elctrical losses described above are already taken into account to get the desired_final_wattage because we account for it when doing the final_cell_amperage

final_cell_amperage * actual_cell_voltage * required_number_of_cells
#above should match desired_final_wattage, which it does! (25773.2 V)

total_waste_heat = stack_chemical_waste_heat + resistive_heat

stack_total_heat_and_electrical_power = desired_final_wattage + total_waste_heat

fuel_conversion_factor = 0.90

fuel_molar_flow_per_cell = uconvert(u"mol/s", final_cell_amperage / (6.0 * fuel_conversion_factor * faraday_constant)) #0.0002313939908035979 mol s^-1
fuel_mass_flow_per_cell = uconvert(u"kg/s", fuel_molar_flow_per_cell * 32.04u"g/mol") #7.413863465347276e-6 kg s^-1

fuel_molar_flow_total = fuel_molar_flow_per_cell * required_number_of_cells #0.07335189508474053 mol s^-1
fuel_mass_flow_total = uconvert(u"kg/s", fuel_molar_flow_total * 32.04u"g/mol") #0.0023501947185150867 kg s^-1



steam_to_carbon_ratio = 1.3
fuel_carbon_ratio = 1

steam_molar_flow_total = steam_to_carbon_ratio * fuel_carbon_ratio * fuel_molar_flow_total

CO_molar_flow_per_cell = 3 * (fuel_molar_flow_per_cell * fuel_conversion_factor)
CO_molar_flow_total = 3 * (fuel_molar_flow_total * fuel_conversion_factor)

H2_molar_flow_per_cell = 7 * (fuel_molar_flow_per_cell * fuel_conversion_factor)
H2_molar_flow_total = 7 * (fuel_molar_flow_total * fuel_conversion_factor)

struct SOFCParameters
    # --- Operating Conditions ---
    temperature::typeof(1.0u"K")
    pressure::typeof(1.0u"Pa")

    # --- Geometric Properties ---
    active_area::typeof(1.0u"m^2") # Area per cell

    # --- Material/Electrochemical Properties ---
    # Ohmic Resistance
    electrolyte_thickness::typeof(1.0u"m")
    electrolyte_conductivity::typeof(1.0u"Ω") # Conductivity as a function of temperature
    ASR_ohmic_other::typeof(1.0u"Ω/m") # Resistance from electrodes, contacts, etc.

    # Activation Losses
    i0_anode::typeof(1.0u"A/m^2")      # Exchange current density for anode
    i0_cathode::typeof(1.0u"A/m^2")    # Exchange current density for cathode
    alpha_anode::Float64              # Anode charge transfer coefficient (typically ~0.5)
    alpha_cathode::Float64            # Cathode charge transfer coefficient (typically ~0.5)
    
    # --- Constants ---
    n_electrons::Int # Number of electrons per H2 molecule
    faraday_const::typeof(1.0u"C/mol")
    gas_const::typeof(1.0u"J/(mol*K)")
end

# --- Example of creating an instance of this struct ---

# Define a function for YSZ electrolyte conductivity (common in SOFCs)
# This is an Arrhenius-type relationship
function ysz_conductivity(T)
    T_K = ustrip(u"K", T)
    sigma_0 = 3.34e4 # S/m
    Ea_cond = 0.83   # eV
    k_boltzmann = 8.617e-5 # eV/K
    return sigma_0 * exp(-Ea_cond / (k_boltzmann * T_K)) * u"S/m"
end

# Now, create a set of parameters for a typical cell
sofc_params = SOFCParameters(
    # Operating Conditions
    750u"°C" |> u"K", #temperature
    1.5u"atm" |> u"Pa", #pressure
    
    # Geometry
    100u"cm^2" |> u"m^2", #active area # e.g., a 10cm x 10cm cell
    
    # Ohmic
    10u"μm" |> u"m", #electrolyte thickness
    100u"Ω",
    #(750u"°C" |> u"K"), #electrolyte conductivity
    0.05u"Ω/cm" |> u"Ω/m", #ASR ohmic factor # ASR for non-electrolyte components
    
    # Activation
    1500u"A/m^2", #anode current exchange density ## These are highly dependent on T and materials
    1000u"A/m^2", #cathode current exchange density # Find values from literature for your materials
    0.5, #anode transfer coefficient
    0.5, #cathode transfer coefficient

    # Constants
    2, #electrons created per mole of fuel reacted (2 for hydrogen and carbon monoxide)
    96485.33u"C/mol",
    8.3145u"J/(mol*K)"
)

function calculate_ohmic_loss(p, j)
    # 1. Calculate ASR of the electrolyte
    conductivity = p.electrolyte_conductivity
    ASR_electrolyte = p.electrolyte_thickness / conductivity # Units: (m) / (S/m) = Ω·m²

    # 2. Total ASR is electrolyte + other components
    ASR_total = ASR_electrolyte + p.ASR_ohmic_other

    # 3. Calculate voltage loss (η = j * ASR)
    eta_ohmic = j * ASR_total
    return uconvert(u"V", eta_ohmic)
end

function calculate_activation_loss(p, j)
    R = p.gas_const
    T = p.temperature
    F = p.faraday_const
    n = p.n_electrons
    
    # Anode Loss
    # We need j > p.i0_anode for this simplification to hold
    eta_act_anode = (R * T / (p.alpha_anode * n * F)) * asinh(j / (2 * p.i0_anode))

    # Cathode Loss
    eta_act_cathode = (R * T / (p.alpha_cathode * n * F)) * asinh(j / (2 * p.i0_cathode))

    total_eta_act = eta_act_anode + eta_act_cathode
    return uconvert(u"V", total_eta_act)
end

function calculate_potentials(p, j, bulk_partial_pressures)
    # bulk_partial_pressures is a Dict: Dict("H2" => 1atm, "H2O" => 0.2atm, "O2" => 0.21atm)
    R = p.gas_const
    T = p.temperature
    F = p.faraday_const
    n = p.n_electrons

    # --- 1. Ideal Nernst Potential (using bulk gas composition) ---
    # Get standard potential (can be a function of T for more accuracy)
    E0 = 1.23u"V" - (T - 298u"K") * 0.00085u"V/K" # Empirical fit for H2/O2
    
    p_h2_bulk = uconvert(u"Pa", bulk_partial_pressures["H2"])
    p_h2o_bulk = uconvert(u"Pa", bulk_partial_pressures["H2O"])
    p_o2_bulk = uconvert(u"Pa", bulk_partial_pressures["O2"])

    # Reference pressure for activity calculation
    p_ref = 1.0u"atm" |> u"Pa"

    E_nernst = E0 - (R * T / (n * F)) * log( (p_h2o_bulk/p_ref) / ((p_h2_bulk/p_ref) * (p_o2_bulk/p_ref)^0.5) )

    # --- 2. Concentration Loss ---
    # This is often modeled using a limiting current density (j_L)
    # j_L is the current at which fuel is consumed faster than it can be supplied
    # This requires a detailed mass transport model. A simpler approach is:
    # Let's assume a limiting current density for now. This should be calculated later.
    jL_anode = 15000u"A/m^2"  # Example value
    jL_cathode = 20000u"A/m^2" # Example value

    eta_conc_anode = (R * T / (n * F)) * log(jL_anode / (jL_anode - j))
    eta_conc_cathode = (R * T / (n * F)) * log(jL_cathode / (jL_cathode - j))
    
    total_eta_conc = eta_conc_anode + eta_conc_cathode

    return uconvert(u"V", E_nernst), uconvert(u"V", total_eta_conc)
end

function solve_cell_voltage(p, j, bulk_partial_pressures)
    # Ensure current is not exceeding the limit
    # (In a real model, j_L would be passed in or calculated inside calculate_potentials)
    if j >= 15000u"A/m^2" 
        return 0.0u"V" # Cell voltage collapses past limiting current
    end
    
    E_nernst, eta_conc = calculate_potentials(p, j, bulk_partial_pressures)
    eta_ohmic = calculate_ohmic_loss(p, j)
    eta_act = calculate_activation_loss(p, j)

    V_cell = E_nernst - eta_act - eta_ohmic - eta_conc
    
    return V_cell
end

# --- Example Usage ---
# Let's test it at a specific current density
target_j = 0.8u"A/cm^2" |> u"A/m^2"

# Assume an inlet gas composition from your reformer
# THIS IS THE KEY LINK TO YOUR OTHER SCRIPT
reformer_outlet_comp = Dict(
    "H2" => 0.5 * sofc_params.pressure,
    "H2O" => 0.15 * sofc_params.pressure,
    "CO" => 0.1 * sofc_params.pressure,
    "CO2" => 0.05 * sofc_params.pressure
)

# Air at the cathode
cathode_gas_comp = Dict("O2" => 0.21 * sofc_params.pressure)

# Combine into one dictionary for the function
partial_pressures = merge(reformer_outlet_comp, cathode_gas_comp)


V_actual = solve_cell_voltage(sofc_params, target_j, partial_pressures)
println("At a current density of $target_j, the predicted cell voltage is $V_actual")

# You can now generate a polarization curve
j_range = (0.01:0.01:1.5)u"A/cm^2"
V_curve = [ustrip(solve_cell_voltage(sofc_params, j |> u"A/m^2", partial_pressures)) for j in j_range]

# using Plots
# plot(ustrip.(j_range), V_curve, xlabel="Current Density (A/cm^2)", ylabel="Cell Voltage (V)", label="Polarization Curve", lw=2)
# title!("SOFC Performance at T=$(sofc_params.temperature)")
