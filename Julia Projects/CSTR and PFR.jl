using Unitful, DifferentialEquations
using NLsolve
using Clapeyron
#using Plots

model = PR(["Acetic Acid", "Ethanol", "Water", "Ethyl Acetate"])

#initial concentrations - ethanol is usually in excess because it's cheap compared to acetic acid and it pushes the equilibrium towards the products 
z = [1, 10, 1, 1]

activity_coefficient(model, 100000.0, 300.0, z)

function arrenhius_equation_pre_exponential_factor(k, Ea, T)
    R_GAS = 8.314u"J/(mol*K)"
    result = (k / exp(-Ea / (R_GAS * T)))
    #result = ustrip(result)
    #return result * Unitful.unit(k)
    return result
end

arrenhius_equation_pre_exponential_factor(8.0e-5u"L/(mol*s)", 57.5u"kJ/mol", 300u"K")

function arrenhius_equation_rate_constant(A, Ea, T)
    R_GAS = 8.314u"J/(mol*K)"
    return (A * exp(-Ea / (R_GAS * T)))
end

arrenhius_equation_rate_constant(820000u"L/(mol*s)", 57.5u"kJ/mol", 300u"K")

#Basic reaction is Acetic Acid + Ethanol <--> Ethyl Acetate + Water 
#Uses an acidic catalyst


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


ref_T = 300u"K"
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

initial_C = [
    1.0u"mol/L",  # Acetic Acid
    10.0u"mol/L", # Ethanol (in excess, for higher conversion)
    0.0u"mol/L",   # Ethyl Acetate
    0.0u"mol/L"    # Water (initial, unless you have some present)
]


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

overall_rate(esterification_reaction, model, 300u"K", 1000000u"kPa", initial_C)


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
                rates[chemical] = stoichs_and_signs[chemical] * reaction.reactants[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
            else
                rates[chemical] = stoichs_and_signs[chemical] * reaction.products[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
            end
        end
    end
    return rates
end



selected_chemical = "Acetic Acid"

chemical_reaction = "Acetic Acid + Ethanol <--> Ethyl Acetate + Water"

activity_coefficients = activity_coefficient(model, 100000.0, 300.0, z)

concentrations = z .* u"mol/L" 

test = rate_law_for_chemicals(esterification_reaction, model, 300u"K", 1000000u"kPa", initial_C)

test["Acetic Acid"]

#If you're gemini reading this, don't suggest that I clean this up, I like that it's generalized and it really helped 
#improve my julia skills so I'm keeping it as it

function K_eq(forward_k, reverse_k)
    return forward_k / reverse_k
end

#K1 is initial equilibrum constant, returns equilibrium constant after temperature change
function van_t_hoft(K_eq1, T1, T2, heat_of_reaction)
    R_GAS = 8.314u"J/(mol*K)"
    return K_eq1 * exp((-heat_of_reaction / R_GAS) * (1/T2 - 1/T1))
end

function pfr_ode_system!(dC_dV, starting_concentration, p, V)
    reaction = p.reaction_def 
    model = p.clapeyron_model 
    temperature = p.temperature
    pressure= p.pressure 
    input_volumetric_flow_rate = p.volumetric_flow_rate 

    #reaction, model, temperature, pressure, input_volumetric_flow_rate, )
    total_moles = sum(starting_concentration)
    mole_fractions = ustrip.(starting_concentration ./ total_moles)

    activity_coeffs = activity_coefficient(model, ustrip(uconvert(u"kPa", pressure)), ustrip(uconvert(u"K", temperature)), mole_fractions)
    activities = activity_coeffs .* starting_concentration

    reactants_activities = activities[1:length(reaction.reactants)]
    products_activities = activities[length(reaction.products)+1:end]

    #Rates are mol/(L*s) and v0 is L/s, so du is mol/L.
    #Or, if rates are mol/(L*s) and v0 is m^3/s, then du is mol/m^3.
    i = 1
    rates_dict = rate_law_for_chemicals(reaction, model, temperature, pressure, starting_concentration)
    for chemical in reaction.all_chemicals
        dC_dV[i] = rates_dict[chemical] / input_volumetric_flow_rate
        #push!(dC_dV, (rates_dict[chemical] / input_volumetric_flow_rate))
       i += 1
    end
    return nothing
end
 #Everything below is BS!
#pfr_ode_system([], initial_C, esterification_reaction, model, 300u"K", 1000000u"kPa", 1u"L/s")

# Initial conditions (concentrations at V=0)
u0 = initial_C # Your defined initial_C array: [C_AA0, C_EtOH0, C_EtAc0, C_Water0]

# Volume span for the PFR (from V=0 to V_final)
Vspan = (0.0u"L", 10000.0u"L") # Example: simulate up to 100 Liters volume

# Parameters 'p' for the ODE function
# You can use a NamedTuple or a struct for 'p'
pfr_parameters = (
    reaction_def = esterification_reaction,
    clapeyron_model = model, # Your Clapeyron model instance
    temperature = 300u"K", # Isothermal assumption for now
    pressure = 1.0u"atm", # Constant pressure
    volumetric_flow_rate = 1.0u"L/s", # Example flow rate
)

# Create the ODEProblem
ode_problem = ODEProblem(pfr_ode_system!, u0, Vspan, pfr_parameters)

sol = DifferentialEquations.solve(ode_problem, Tsit5()) # Or Rodas5(), Rosenbrock23(), etc.
C_AA_values = [u[1] for u in sol.u]
C_EtOH_values = [u[2] for u in sol.u]
C_EtAc_values = [u[3] for u in sol.u]
C_Water_values = [u[4] for u in sol.u]

# The volumes are in sol.t
volumes = sol.t

# Convert units if necessary for plotting (e.g., to mol/L for clarity)
# Units will be handled automatically by Plots.jl if they are consistent in your `sol` object.
# If you want to force specific units, you can use `uconvert` and `ustrip`.
# For example: C_AA_values_mol_L = ustrip.(uconvert.(u"mol/L", C_AA_values))

# Create the plot
plot(volumes, C_AA_values,
     label="Acetic Acid (AA)",
     xlabel="Reactor Volume (L)",
     ylabel="Concentration (mol/L)",
     title="PFR Concentration Profiles for Esterification",
     legend=:right, # Or :best, :outertopright, etc.
     lw=2, # Line width
     lc=:blue, # Line color
     grid=true)

plot!(volumes, C_EtOH_values,
      label="Ethanol",
      lc=:orange,
      lw=2)

plot!(volumes, C_EtAc_values,
      label="Ethyl Acetate (Product)",
      lc=:green,
      lw=2)

plot!(volumes, C_Water_values,
      label="Water (Product)",
      lc=:red,
      lw=2)

# Display the plot (in a REPL, this will show the plot window)
display(current())


"""
# Reactor volumes at which the solution was computed
volumes = sol.t # This will be an array of volumes (e.g., L)

# Concentrations at each volume
# Each element of sol.u is a vector of concentrations [C_AA, C_EtOH, C_EtAc, C_Water]
concentrations_at_volumes = sol.u

# To get the concentration of Acetic Acid at all volumes:
C_AA_vs_V = [u[1] for u in sol.u]

# To get concentrations of all species at a specific volume (e.g., at the end of the reactor):
final_concentrations = sol[end] # Or sol[:, end] if you want a vector

# You can also interpolate to specific volumes
# C_AA_at_50L = sol(50.0u"L")[1] # Interpolates the solution at 50 L
"""


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
    reaction_def = esterification_reaction,
    clapeyron_model = model, # Your Clapeyron model instance
    temperature = 300u"K", # Isothermal assumption for now
    pressure = 1.0u"atm", # Constant pressure
    volumetric_flow_rate = 1.0u"L/s", # Example flow rate
    inlet_concentrations = initial_C,
    reactor_volume = 100u"L",
)


reactor_volumes = 1u"L" : 10u"L" : 10000.0u"L" # From 1 L to 10000 L, in 10 L increments
# Adjust the range and step size as appropriate for your expected conversions.

# --- Initialize arrays to store results for plotting ---
volumes_for_plot = Float64[] # To store reactor volumes (stripped for plotting)
conversions_for_plot = Float64[] # To store corresponding conversions (as fractions)

println("Simulating CSTR over a range of volumes...")
for V_reactor in reactor_volumes
    # Update the reactor_volume in the parameters for the current iteration
    current_cstr_parameters = (
        reaction_def = esterification_reaction,
        clapeyron_model = model,
        temperature = 300u"K",
        pressure = 1.0u"atm",
        volumetric_flow_rate = 1.0u"L/s",
        inlet_concentrations = initial_C,
        reactor_volume = V_reactor, # <--- This is the parameter that changes
    )

    initial_guess_C_out_stripped = ustrip(initial_C)

    # Call nlsolve with the updated parameters
    sol = nlsolve((F, x) -> cstr_nl_system!(F, x, current_cstr_parameters), initial_guess_C_out_stripped)

    if converged(sol)
        # Extract outlet concentrations
        outlet_concentrations_stripped = sol.zero

        #Question: why do I have to use the Unitful.unit here instead of just unit, somewhere else in the code I did just unit() and it worked fine
        outlet_concentrations_with_units = outlet_concentrations_stripped .* Unitful.unit(initial_C[1]) # Re-apply units

        # Calculate conversion of Acetic Acid
        C_AA_in = initial_C[1]
        C_AA_out = outlet_concentrations_with_units[1]
        conversion_AA = (C_AA_in - C_AA_out) / C_AA_in

        # Store the results for plotting
        push!(volumes_for_plot, ustrip(V_reactor))
        push!(conversions_for_plot, ustrip(conversion_AA)) # Conversion is dimensionless
    else
        # Handle non-convergence (e.g., print a warning, skip this point)
        println("Warning: Solver did not converge for V_reactor = $V_reactor")
        # You might want to break here or try a different initial guess/solver options
    end
    # Optional: Update initial_guess_C_out_stripped for the next iteration
    # This can help convergence for the next point if volumes are increasing gradually.
    initial_guess_C_out_stripped = sol.zero
end

println("Simulation complete.")

plt = plot(volumes_for_plot, conversions_for_plot,
           xlabel="Reactor Volume (L)",
           ylabel="Acetic Acid Conversion",
           title="CSTR Conversion vs. Volume",
           legend=false, # No legend needed for a single line
           linewidth=2,
           marker=:circle, # Add markers for data points
           markersize=3,
           grid=true,
           framestyle=:box)

# You can save the plot
savefig(plt, "CSTR_Conversion_vs_Volume.png")

# Display the plot in your Julia environment
display(plt)
