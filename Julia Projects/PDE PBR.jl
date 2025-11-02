using Revise

using Unitful, DifferentialEquations
using Clapeyron
using Plots

using ModelingToolkit
using MethodOfLines
using DomainSets

function arrenhius_equation_pre_exponential_factor(k, Ea, T)
    R = 8.314u"J/(mol*K)"
    result = (k / exp(-Ea / (R * T)))
    #result = ustrip(result)
    #return result * Unitful.unit(k)
    return result
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

function van_t_hoft(K_eq1, T1, T2, heat_of_reaction)
    R = 8.314u"J/(mol*K)"
    return K_eq1 * exp((-heat_of_reaction / R) * (1/T2 - 1/T1))
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

reactor_inlet_temperature = 300.13u"°C" |> u"K"

reactor_inlet_pressure = 1.0u"bar"

SMR_methanol_molar_flow = 0.001u"mol/s"

SMR_model = PR(["Methanol", "water", "carbon monoxide", "hydrogen", "Carbon Dioxide"])

mole_fractions = [1.0, 1.3, 0.001, 0.001, 0.001]
SMR_initial_C = mole_fractions .* u"mol/L"

molar_flow_correction_factor = SMR_methanol_molar_flow / SMR_initial_C[1]

SMR_molar_flows = SMR_initial_C .* molar_flow_correction_factor

total_molar_flow = sum(SMR_molar_flows)

# CH3OH + H2O ⇋ CO + 3H2
# Below is for SMR
ref_T = reactor_inlet_temperature
kf_ref = 1e-1u"mol/(kg*s*bar^2)" 
Ea_f = 90u"kJ/mol" 
kf_A = arrenhius_equation_pre_exponential_factor(kf_ref, Ea_f, ref_T)

Keq_ref = 100000.0
#K_gibbs_free(298u"K", reactor_inlet_temperature, 33.2u"kJ/mol", 49.4u"kJ/mol")
#formatted reference temperature, actual temperature, delta gibbs free energy of reaction, and heat of reaction

kr_ref = kf_ref / Keq_ref
# For Ea_r, you know Ea_r = Ea_f - dH_rxn, assuming elementary.
Ea_r = 43u"kJ/mol"
kr_A = arrenhius_equation_pre_exponential_factor(kr_ref, Ea_r, ref_T) / 1.0u"bar^2" 
#The Above /1.0u"bar^s" is very necessary because it causes the net_eraction_rate to not error on a dimension error (trust me!)

kf_A_SMR_val = ustrip(kf_A)
kr_A_SMR_val = ustrip(kr_A)
Ea_f_SMR_val = ustrip(Ea_f)
Ea_r_SMR_val = ustrip(Ea_r)
Keq_ref_SMR_val = ustrip(Keq_ref)


SMR_reaction = HeterogeneousChemicalReaction(
    "SMR",
    Dict("CH3OH" => 1, "H2O" => 1),
    Dict("CO" => 1, "H2" => 3),
    ["CH3OH", "H2O", "CO", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)

#rate_law_for_chemicals(SMR_reaction, SMR_model, reactor_inlet_temperature, 100000u"Pa", SMR_initial_C)

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

kf_A_WGS_val = ustrip(kf_A)
kr_A_WGS_val = ustrip(kr_A)
Ea_f_WGS_val = ustrip(Ea_f)
Ea_r_WGS_val = ustrip(Ea_r)
Keq_ref_WGS_val = ustrip(Keq_ref)

F_CH3OH_in = ustrip(0.001u"mol/s")
F_H2O_in = F_CH3OH_in * 1.3
F_CO_in = F_CH3OH_in * 0.001 # Assuming small initial amount
F_H2_in = F_CH3OH_in * 0.001
F_CO2_in = F_CH3OH_in * 0.001

MW_CH3OH_val = 0.03204
MW_H2O_val  = 0.018015
MW_CO_val   = 0.02801
MW_H2_val   = 0.002016
MW_CO2_val  = 0.04401

@parameters t z 

@variables R_GAS A_c ρ_bulk Dp epsilon mu_gas MW_CH3OH MW_H2O MW_CO MW_H2 MW_CO2 kf_A_SMR Ea_f_SMR kr_A_SMR Ea_r_SMR Keq_ref_SMR kf_A_WGS Ea_f_WGS kr_A_WGS Ea_r_WGS Keq_ref_WGS

@syms F_CH3OH(t, z)::Real F_H2O(t, z)::Real F_CO(t, z)::Real F_H2(t, z)::Real F_CO2(t, z)::Real T(t, z)::Real P(t, z)::Real

Dt = Differential(t)
Dz = Differential(z)

F_total = F_CH3OH(t, z) + F_H2O(t, z) + F_CO(t, z) + F_H2(t, z) + F_CO2(t, z)

y_CH3OH = F_CH3OH(t, z) / F_total 
y_H2O = F_H2O(t, z) / F_total 
y_CO = F_CO(t, z) / F_total 
y_H2 = F_H2(t, z) / F_total 
y_CO2 = F_CO2(t, z) / F_total

p_CH3OH = y_CH3OH * P(t, z)
p_H2O = y_H2O * P(t, z)
p_CO = y_CO * P(t, z)
p_H2 = y_H2 * P(t, z)
p_CO2 = y_CO2 * P(t, z)

kf_SMR = kf_A_SMR * exp(-Ea_f_SMR / (R_GAS * T(t, z)))
kr_SMR = kf_SMR / Keq_ref_SMR

kf_WGS = kf_A_WGS * exp(-Ea_f_WGS / (R_GAS * T(t, z)))
kr_WGS = kf_WGS / Keq_ref_WGS

r_net_SMR = kf_SMR * (p_CH3OH * p_H2O) - kr_SMR * (p_CO * p_H2^3)
r_net_WGS = kf_WGS * (p_CO * p_H2O) - kr_WGS * (p_CO2 * p_H2)

eq_F_CH3OH = Dz(F_CH3OH(t, z)) ~ ρ_bulk * A_c * (-1 * r_net_SMR)
eq_F_H2O = Dz(F_H2O(t, z)) ~ ρ_bulk * A_c * (-1 * r_net_SMR + -1 * r_net_WGS)
eq_F_CO = Dz(F_CO(t, z)) ~ ρ_bulk * A_c * (1 * r_net_SMR + -1 * r_net_WGS)
eq_F_H2 = Dz(F_H2(t, z)) ~ ρ_bulk * A_c * (3 * r_net_SMR + 1 * r_net_WGS)
eq_F_CO2 = Dz(F_CO2(t, z)) ~ ρ_bulk * A_c * (1 * r_net_WGS)

MW_avg = y_CH3OH * MW_CH3OH + y_H2O * MW_H2O + y_CO * MW_CO + y_H2 * MW_H2 + y_CO2 * MW_CO2
mass_flow_total = F_CH3OH(t, z) * MW_CH3OH + F_H2O(t, z) * MW_H2O + F_CO(t, z) * MW_CO + F_H2(t, z) * MW_H2 + F_CO2(t, z) * MW_CO2
G = mass_flow_total / A_c
rho_gas = (P(t, z) * MW_avg) / (R_GAS * T(t, z))
ergun_term_1 = 150 * mu_gas * (1 - epsilon)^2 / (Dp^2 * epsilon^3)
ergun_term_2 = 1.75 * G * (1 - epsilon) / (Dp * epsilon^3)
ergun_expression = - (G / rho_gas) * (ergun_term_1 + ergun_term_2)

# Your pressure equation is then:
eq_P = Dz(P(t, z)) ~ ergun_expression

eq_T = Dz(T(t, z)) ~ 0

eqs = [eq_F_CH3OH, eq_F_H2O, eq_F_CO, eq_F_H2, eq_F_CO2, eq_P, eq_T]

domains = [t ∈ Interval(0.0, 100.0),
           z ∈ Interval(0.0, 1.0)]  

test(t, z) = 0.001 / (100 - 0)

bcs = [
    F_CH3OH(0, z) ~ 0.0,
    F_H2O(t, 0) ~ F_H2O_in,
    F_CO(t, 0) ~ F_CO_in,
    F_H2(t, 0) ~ F_H2_in,
    F_CO2(t, 0) ~ F_CO2_in,
    P(t, 0) ~ ustrip(reactor_inlet_pressure),
    T(t, 0) ~ ustrip(reactor_inlet_temperature)
]

@named pdesystem = PDESystem(eqs, bcs, domains, [t, z], [F_CH3OH(t, z), F_H2O(t, z), F_CO(t, z), F_H2(t, z), F_CO2(t, z), T(t, z), P(t, z)])

pdesystem.domain
pdesystem.bcs


discretization = MOLFiniteDifference([z => 20], t) # Discretize z into 20 points
prob = discretize(pdesystem, discretization)

# Solve the problem
# Since this is a steady-state problem, the time `t` is just a dummy variable.
# The solver integrates from z=0 to z=reactor_length.
sol = solve(prob, Tsit5())

# --- Plotting Results ---
# We need to extract the solution for each variable at different points in z
z_grid = sol.ivs[2]
F_CH3OH_sol = sol[F_CH3OH(t, z)]
P_sol = sol[P(t, z)]
# etc. for other variables

# Plot pressure drop
plot(z_grid, P_sol ./ 1e5, xlabel="Reactor Length (m)", ylabel="Pressure (bar)", label="Pressure", lw=2)

# Plot molar flows
plot(z_grid, sol[F_CH3OH(t, z)], xlabel="Reactor Length (m)", ylabel="Molar Flow (mol/s)", label="CH3OH", lw=2)
plot!(z_grid, sol[F_H2(t, z)], label="H2")
plot!(z_grid, sol[F_CO(t, z)], label="CO")