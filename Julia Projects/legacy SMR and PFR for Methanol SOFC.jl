
using Revise

using Unitful, DifferentialEquations
using NLsolve
using Clapeyron
using Plots
using Symbolics
#include("Chemeqs/src/Chemeqs.jl")
using .Chemeqs 
#import .Chemeqs:overall_rate, arrenhius_equation_pre_exponential_factor, rate_law_for_chemicals #For some reason this only works if I do this

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
SMR_initial_mole_fractions = z .* u"mol"

# CH3OH + H2O ⇋ CO + 3H2
#K_gibbs_free(temporary_reaction, SMR_model, reactor_temperature, 1u"atm", z, 49.4u"kJ/mol")
# Below is for SMR
ref_T = reactor_temperature
kf_ref = 8.0e-5u"L/(mol*s)"
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

SMR_reaction = ChemicalReaction(
    "SMR",
    Dict("CH3OH" => 1, "H2O" => 1),
    Dict("CO" => 1, "H2" => 3),
    ["CH3OH", "H2O", "CO", "H2"],
    kf_A, Ea_f,
    kr_A, Ea_r
)


#overall_rate(SMR_reaction, SMR_model, reactor_temperature, 100000u"Pa", SMR_initial_C)

rate_law_for_chemicals(SMR_reaction, SMR_model, reactor_temperature, 100000u"Pa", SMR_initial_C)







function pfr_ode_system!(dC_dV, C, p, V)
    reaction = p.reaction_def 
    model = p.clapeyron_model 
    temperature = p.temperature
    pressure= p.pressure 
    input_volumetric_flow_rate = p.volumetric_flow_rate 

    #reaction, model, temperature, pressure, input_volumetric_flow_rate, )
    total_moles = sum(C)
    mole_fractions = ustrip.(C ./ total_moles)

    activity_coeffs = activity_coefficient(model, ustrip(uconvert(u"Pa", pressure)), ustrip(uconvert(u"K", temperature)), mole_fractions)
    activities = activity_coeffs .* C

    reactants_activities = activities[1:length(reaction.reactants)]
    products_activities = activities[length(reaction.products)+1:end]

    #Rates are mol/(L*s) and v0 is L/s, so du is mol/L.
    #Or, if rates are mol/(L*s) and v0 is m^3/s, then du is mol/m^3.
    i = 1
    rates_dict = rate_law_for_chemicals(reaction, model, temperature, pressure, C)
    for chemical in reaction.all_chemicals
        dC_dV[i] = rates_dict[chemical] / input_volumetric_flow_rate
        #push!(dC_dV, (rates_dict[chemical] / input_volumetric_flow_rate))
       i += 1
    end
    return nothing
end

# Initial conditions (concentrations at V=0)
u0 = SMR_initial_C # Your defined initial_C array: [C_Methanol0, C_H2O0, C_CO0, C_H20]

# Volume span for the PFR (from V=0 to V_final)
Vspan = (1.0u"L", 600000.0u"L") # For some reason I only get good conversion if i make the total reactor volume that huge

# Parameters 'p' for the ODE function
# You can use a NamedTuple or a struct for 'p'
methanol_model = PR(["Methanol"])
pfr_parameters = (
    reaction_def = SMR_reaction,
    clapeyron_model = SMR_model, # Your Clapeyron model instance
    temperature = reactor_temperature, # Isothermal assumption for now
    pressure = 1.0u"atm", # Constant pressure
    volumetric_flow_rate = SMR_methanol_molar_flow * volume(methanol_model, 100000, ustrip(reactor_temperature), [1.0]) * u"m^3/mol"
    #1.0u"L/s"
    #8.544543774765826e-5u"L/s", 
)

# Create the ODEProblem
ode_problem = ODEProblem(pfr_ode_system!, u0, Vspan, pfr_parameters)

pfr_sol = DifferentialEquations.solve(ode_problem, Tsit5()) # Or Rodas5(), Rosenbrock23(), etc.
C_Methanol_values = [u[1] for u in pfr_sol.u]
C_H2O_values = [u[2] for u in pfr_sol.u]
C_CO_values = [u[3] for u in pfr_sol.u]
C_H2_values = [u[4] for u in pfr_sol.u]

heat_required = uconvert(u"W", ((C_Methanol_values[1] - C_Methanol_values[end]) / C_Methanol_values[1]) * 1 * SMR_methanol_molar_flow * 49.4u"kJ/mol")

# The volumes are in pfr_sol.t
volumes = pfr_sol.t

plt = plot(volumes, C_Methanol_values,
    xlabel="Reactor Volume",
    ylabel="Methanol Concentration",
    title="PFR Methnaol Concentration vs.\nVolume and Heat Input",
    legend=:topleft, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    label="Methanol Concentration mol/L",
    framestyle=:box)

wattages_required = []
for i in range(1, length(C_Methanol_values))
    push!(wattages_required, uconvert(u"W", ((C_Methanol_values[1] - C_Methanol_values[i]) / C_Methanol_values[1]) * 1 * SMR_methanol_molar_flow * 49.4u"kJ/mol") )
end
wattages_required = ustrip.(wattages_required)

plt = plot!(twinx(plt), volumes, wattages_required,
    ylabel="Heat Input (W)",
    legend=:bottomleft, # No legend needed for a single line
    linewidth=2,
    marker=:circle, # Add markers for data points
    markersize=3,
    grid=true,
    color=:red,
    label="Heat Input (W)",
    framestyle=:box)








function cstr_nl_system!(F, outlet_concentations, p)
    reaction = p.reaction_def 
    model = p.clapeyron_model 
    temperature = p.temperature
    pressure= p.pressure 
    inlet_concentrations = p.inlet_concentrations 
    volumetric_flow_rate = p.volumetric_flow_rate
    reactor_volume = p.reactor_volume 

    rates_dict = rate_law_for_chemicals(reaction, model, temperature, pressure, outlet_concentations)
    i = 1
    for chemical in reaction.all_chemicals
        molar_flow_of_chem_in = ustrip(inlet_concentrations[i] * volumetric_flow_rate)
        molar_flow_of_chem_out = ustrip(outlet_concentations[i] * volumetric_flow_rate)
        F[i] = (molar_flow_of_chem_in - molar_flow_of_chem_out) + (ustrip(rates_dict[chemical]) * ustrip(reactor_volume))
        i += 1
    end
    return nothing
end

cstr_parameters = (
    reaction_def = SMR_reaction,
    clapeyron_model = SMR_model, # Your Clapeyron model instance
    temperature = reactor_temperature, # Isothermal assumption for now
    pressure = 1.0u"atm", # Constant pressure
    volumetric_flow_rate = SMR_methanol_molar_flow * volume(methanol_model, 100000, ustrip(reactor_temperature), [1.0]) * u"m^3/mol",
    #1.0u"L/s", # Example flow rate
    inlet_concentrations = SMR_initial_C,
    reactor_volume = 100.0u"L",
)


reactor_volumes = 1.0u"L" : 50u"L" : 2000.0u"L" #Formatted start, step, stop
# Adjust the range and step size as appropriate for your expected conversions.

# --- Initialize arrays to store results for plotting ---
volumes_for_plot = Float64[] # To store reactor volumes (stripped for plotting)
conversions_for_plot = Float64[] # To store corresponding conversions (as fractions)

println("Simulating CSTR over a range of volumes...")
for V_reactor in reactor_volumes
    # Update the reactor_volume in the parameters for the current iteration
    current_cstr_parameters = (
        reaction_def = SMR_reaction,
        clapeyron_model = SMR_model, # Your Clapeyron model instance
        temperature = reactor_temperature, # Isothermal assumption for now
        pressure = 1.0u"atm", # Constant pressure
        volumetric_flow_rate = SMR_methanol_molar_flow * volume(methanol_model, 100000, ustrip(reactor_temperature), [1.0]) * u"m^3/mol",
        inlet_concentrations = SMR_initial_C,
        reactor_volume = V_reactor,
    )

    initial_guess_C_out_stripped = ustrip(SMR_initial_C)

    # Call nlsolve with the updated parameters
    cstr_sol = nlsolve((F, x) -> cstr_nl_system!(F, x, current_cstr_parameters), initial_guess_C_out_stripped)

    if converged(cstr_sol)
        # Extract outlet concentrations
        outlet_concentrations_stripped = cstr_sol.zero

        #Keep the Unitful.unit
        outlet_concentrations_with_units = outlet_concentrations_stripped .* Unitful.unit(SMR_initial_C[1]) # Re-apply units

        # Calculate conversion of Acetic Acid
        C_Methanol_in = SMR_initial_C[1]
        C_Methanol_out = outlet_concentrations_with_units[1]
        conversion_AA = (C_Methanol_in - C_Methanol_out) / C_Methanol_in

        # Store the results for plotting
        push!(volumes_for_plot, ustrip(V_reactor))
        push!(conversions_for_plot, ustrip(conversion_AA)) # Conversion is dimensionless
    else
        # Handle non-convergence (e.g., print a warning, skip this point)
        println("Warning: Solver did not converge for V_reactor")
        # You might want to break here or try a different initial guess/solver options
    end
    # Optional: Update initial_guess_C_out_stripped for the next iteration
    # This can help convergence for the next point if volumes are increasing gradually.
    initial_guess_C_out_stripped = cstr_sol.zero
end


heat_required_cstr = uconvert(u"W", ((C_Methanol_values[1] - C_Methanol_values[end]) / C_Methanol_values[1]) * 1 * SMR_methanol_molar_flow * 49.4u"kJ/mol")


println("Simulation complete.")
conversions_for_plot
plt = plot(volumes_for_plot, conversions_for_plot,
           xlabel="Reactor Volume (L)",
           ylabel="Methanol Conversion",
           title="CSTR Conversion vs. Volume and Heat Required",
           legend=:false, # No legend needed for a single line
           linewidth=2,
           marker=:circle, # Add markers for data points
           markersize=3,
           grid=true,
           label="Methanol Conversion",
           framestyle=:box)

wattages_for_plot = conversions_for_plot .* SMR_methanol_molar_flow * 49.4u"kJ/mol"

wattages_for_plot = uconvert.(u"W", wattages_for_plot)

plt = plot!(twinx(plt), volumes_for_plot, wattages_for_plot,
           ylabel="Heat Input",
           legend=:bottomright, # No legend needed for a single line
           linewidth=2,
           marker=:circle, # Add markers for data points
           markersize=3,
           grid=true,
           label="Methanol Conversion and Heat Input",
           framestyle=:box)
