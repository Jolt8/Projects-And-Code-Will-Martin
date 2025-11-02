
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


