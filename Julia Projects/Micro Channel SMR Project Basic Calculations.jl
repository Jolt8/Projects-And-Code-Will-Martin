using Unitful
using CoolProp

# Constants
S_C_ratio = 3.0
methane_LHV = 802.0u"kJ/mol" # Use LHV for high-temp processes
furnace_efficiency = 0.95 # Realistic estimate (flue gas losses)

T_in = 21.13u"°C" |> u"K"
T_out = 700.0u"°C" |> u"K"
P = 1.0u"bar" # 1 bar in Pa

# 1. Heat required for Water (Inlet liquid -> Outlet 700C gas)
h_w_in = PropsSI("HMOLAR", "T", T_in, "P", P, "Water")
h_w_out = PropsSI("HMOLAR", "T", T_out, "P", P, "Water")
ΔH_water = (h_w_out - h_w_in) |> u"kJ/mol"

# 2. Heat required to preheat Methane (Inlet -> Outlet 700C)
h_m_in = PropsSI("HMOLAR", "T", T_in, "P", P, "Methane")
h_m_out = PropsSI("HMOLAR", "T", T_out, "P", P, "Methane")
ΔH_methane_preheat = (h_m_out - h_m_in) |> u"kJ/mol"

# 3. Reaction Energy
# CH4 + H2O -> CO + 3H2  (+206 kJ/mol)
# CO + H2O -> CO2 + H2   (-41 kJ/mol)
ΔH_reforming = 206.0u"kJ/mol" 
ΔH_wgs = -41.0u"kJ/mol"
fraction_wgs = 0.80

total_reaction_energy = ΔH_reforming + (fraction_wgs * ΔH_wgs)

# 4. Total energy required to process 1 mol of Methane
# (1 mol CH4 preheat + 3 mol Water heat + Reaction energy)
total_heat_needed_per_mol_reformed = ΔH_methane_preheat + (S_C_ratio * ΔH_water) + total_reaction_energy

# 5. Energy available from 1 mol of burned methane
heat_available_per_mol_burned = methane_LHV * furnace_efficiency

# 6. Final Ratio
ratio = heat_available_per_mol_burned / total_heat_needed_per_mol_reformed

println("Heat for Steam (per mol CH4): ", S_C_ratio * ΔH_water)
println("Heat for Reaction + Preheat: ", total_reaction_energy + ΔH_methane_preheat)
println("Moles reformed per mole burned: ", ratio)