using Unitful
using Symbolics, Nemo, Groebner
using Clapeyron
using NLsolve

"""
model = Wilson(["Methane", "Hydrogen", "Carbon Monoxide", "Carbon Dioxide", "Water", "Carbon"])

z = [
    1.0u"mol/L",
    10.0u"mol/L",
    0.0u"mol/L",
    0.0u"mol/L",
    0.0u"mol/L"
]

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


ref_T = 950u"°C"
kf_ref = 8.0e-5u"L/(mol*s)"
Ea_f = 57.5u"kJ/mol"
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

# For reverse reaction, calculate kr from Keq. Let's assume Keq = 4 at 300K
Keq_ref = 4.0u"NoUnits" # Keq is dimensionless
kr_ref = kf_ref / Keq_ref
# For Ea_r, you know Ea_r = Ea_f - dH_rxn, assuming elementary.
# dH_rxn = +33.8 kJ/mol (from our previous discussion)
Ea_r = Ea_f - 33.8u"kJ/mol" # Be careful with units here! Ea_f and dH_rxn need to be in same units.
kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T)


esterification_reaction = ChemicalReaction(
    "Esterification",
    Dict("Acetic Acid" => 1, "Ethanol" => 1),
    Dict("Water" => 1, "Ethyl Acetate" => 1),
    ["Acetic Acid", "Ethanol", "Water", "Ethyl Acetate"],
    kf_A, Ea_f,
    kr_A, Ea_r
)

z = [
    1.0u"mol/L",  # Acetic Acid
    10.0u"mol/L", # Ethanol (in excess, for higher conversion)
    0.0u"mol/L",   # Ethyl Acetate
    0.0u"mol/L"    # Water (initial, unless you have some present)
]
"""

#Carbon_balance should be 1

@variables N_H2, N_CO, N_CO2, N_CH4, N_H2O_prod
@variables S_to_C_ratio, O_to_C_ratio, Carbon_balance
@variables K_WGS

mole_balances = [
    N_H2 ~ 2 * N_CO,
    Carbon_balance ~ N_CO + N_CO2,
    4 + 2 * S_to_C_ratio ~ 2 * N_H2 + 2 * N_H2O_prod, #4 * N_CH4 + 2 * S_to_C_ratio ~ N_CO + 2*N_CO2 + N_H2O_prod,
    2 * O_to_C_ratio + S_to_C_ratio ~ N_CO + 2 * N_CO2 + N_H2O_prod,
    K_WGS ~ (N_CO2 * N_H2) / (N_CO * N_H2O_prod)
]


known_values = Dict(
    K_WGS          => 0.364,
    S_to_C_ratio   => 0.8
)

# 4. Variables to solve for
#    We want to find N_H2, N_CO, N_CO2, N_H2O_prod, and O_to_C_ratio
variables_to_solve = [N_H2, N_CO, N_CO2, N_H2O_prod, O_to_C_ratio]

# Step 1: Substitute known values
eqs_subs = Symbolics.substitute.(mole_balances, (known_values,))

# Step 2: Solve N_H2 and N_CO2 in terms of N_CO
sol_N_H2 = solve_for(eqs_subs[1], N_H2) # N_H2 ~ 2*N_CO
sol_N_CO2 = solve_for(eqs_subs[2], N_CO2) # N_CO2 ~ 1 - N_CO

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
sol_N_CO = solve_for(h_balance_simplified, N_CO)
# Now we have N_CO. Let's evaluate it numerically
N_CO_val = Symbolics.value(sol_N_CO[])

# Calculate the rest of the values numerically
N_H2_val = 2 * N_CO_val
N_CO2_val = 1 - N_CO_val
N_H2O_prod_val = 2 * (1 - N_CO_val) / Symbolics.value(known_values[K_WGS])

# Finally, substitute all knowns into the O Balance to solve for O_to_C_ratio
o_balance_simplified = Symbolics.substitute(eqs_subs[4], Dict(N_CO => N_CO_val, N_CO2 => N_CO2_val, N_H2O_prod => N_H2O_prod_val))
sol_O_to_C_ratio = solve_for(o_balance_simplified, O_to_C_ratio)
O_to_C_ratio_val = Symbolics.value(sol_O_to_C_ratio[])

# Print results
println("--- Calculated Molar Flows (per 1 mole CH4 feed) ---")
println("N_CO: ", N_CO_val)
println("N_H2: ", N_H2_val)
println("N_CO2: ", N_CO2_val)
println("N_H2O_prod: ", N_H2O_prod_val)
println("-----------------------------------------------------")
println("S/C Ratio (input): ", Symbolics.value(known_values[S_to_C_ratio]))
println("O/C Ratio (output): ", O_to_C_ratio_val)




#equations below are related to jet velocity from orifice equation
#jet velocity = sqrt((2 * ΔP) / fluid_density)

function jet_velocity_from_orifice(pressure_drop, fluid_density)
    return sqrt((2 * delta_p) / fluid_density)
end

function pressure_drop_from_orifice(orifice_jet_velocity, fluid_density)
    return (orifice_jet_velocity^2 * fluid_density) / 2
end

function required_density_for_pressure_drop_and_velocity(orifice_jet_velocity, pressure_drop)
    (orifice_jet_velocity^2)/(2 * pressure_drop)
end

#Mixing time equations
#time_mixing ~ nozzle diameter / jet velocity 
function mixing_time(nozzle_diameter, jet_velocity)
    return nozzle_diameter / jet_velocity
end

function nozzle_diameter_for_required_mixing_time_and_jet_velocity(mixing_time, jet_velocity)
    return mixing_time * jet_velocity
end

function jet_velocity_for_required_mixing_time_and_nozzle_diameter(mixing_time, nozzle_diameter)
    return mixing_time / nozzle_diameter
end

#Chemical (ignition) time
#time reacting ≈ 1 / (arrhenius rate constant of dominant flame reaction * concentration of reactants^order)



#RULE OF THUMB: time to mix usually has to be lower to be time to react or else the flame lifts off
"""
    - Energy Balance (Axial direction assuming no radial gradients)
    - z here is the coordinate of one of the variables below along the reactor length
    - For example, temperature could be high at z = 1, but lower at z = 10

    - Because it's derivitive, it tells you how the temperature evolves due to the net effect of exothermic and endothermic reactions
    - For this one, you could have two reactions with different reaction rates, such as an exothermic delta_H with high rate of reaction for POX 
        and an endothermic reaction with low rate of reaction for steam methane reforming
"""

@variables z

@syms T(z) m(z) Cp(z)

@syms r1(z) r2(z) delta_h1(z) delta_h2(z)

dTdz = Differential(z)(T(z))
energy_balance_eqn = dTdz ~ -(r1(z)*delta_h1(z) + r2(z)*delta_h2(z)) / (m(z) * Cp(z))

solve_for(energy_balance_eqn, r1(z))



#Next equation set is for choked flow conditions (for gas nozzle at sonic velocity)

#Eqn: mass flow = discharge coefficient * nozzle throat area * stagnation (inlet) pressure * sqrt(specific_heat_ratio / (R_gas * stagnation (inlet) temperature)) * (2 / (specific_heat_ratio + 1))^((specific_heat_ratio + 1))/2(specific_heat_ratio - 1)
@variables mass_flow discharge_coefficient nozzle_throat_area inlet_pressure inlet_temperature Cp Cv
#note that Cp (specific heat at constant pressure) / Cv (specific heat at constant volume) = specific heat ratio
#also, discharge_coefficient is dimensionless and usually around 0.9-1.0

R_gas = 8.314
#u"J/(mol*K)"

choked_flow_eqn = mass_flow ~ discharge_coefficient * nozzle_throat_area * inlet_pressure * sqrt((Cp / Cv) / (R_gas * inlet_temperature)) * (2 / ((Cp / Cv) + 1))^( (((Cp / Cv) + 1)) / (2*((Cp / Cv) - 1)) )

symbolic_solve(choked_flow_eqn, nozzle_throat_area)[1]
