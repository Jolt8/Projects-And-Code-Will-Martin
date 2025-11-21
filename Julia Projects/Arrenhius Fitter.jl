using GLM, DataFrames
struct ChemicalReaction
    name::String
    reactants::Dict{String, Int}
    products::Dict{String, Int}
    all_chemicals::Vector{String}
    kf_A # Pre-exponential factor for forward reaction
    kf_Ea # Activation energy for forward reaction
    kr_A # Pre-exponential factor for reverse reaction
    kr_Ea  # Activation energy for reverse reaction
end

reaction = ChemicalReaction(
    "Ethanol and Acetic Acid Esterification",
    Dict("CH3CH2OH" => 1, "CH3COOH" => 1),
    Dict("CH3COOCH2CH3" => 1, "H2O" => 1),
    ["CH3CH2OH", "CH3COOH", "CH3COOCH2CH3", "H2O"],
    nothing, nothing,
    nothing, nothing
)

reaction_temperature = 45.13u"°C" |> u"K"

species_starting_concentration = 1.0u"mol/L"

time_concentration_pairs = [[0.0u"s", 1.0u"mol/L"], [1.0u"s", 0.5u"mol/L"]]

species_i_for_pairs = 1

reaction_order = 0

for (reactant, nu) in reaction.reactants
    reaction_order += nu
end

reaction_order #this is kinda bad because reactions with more than 1 step can be a different order than it appears based on its stoichiometry

rate_constants = []

reaction_order = 0
trials = Dict(
    200.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 1.884u"mol/L"]],
    300.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 1.679u"mol/L"]],
    350.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 1.379u"mol/L"]]
)

#=
reaction_order = 1
trials = Dict(
    250.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 2.0u"mol/L"]],
    300.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 1.083u"mol/L"]],
    350.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 0.298u"mol/L"]]
)
=#
#=
reaction_order = 2
trials = Dict(
    250.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 1.0u"mol/L"]],
    300.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 0.739u"mol/L"]],
    350.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 0.383u"mol/L"]]
)
=#
#=
reaction_order = 2
trials = Dict(
    250.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 1.0u"mol/L"]],
    300.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 0.739u"mol/L"]],
    350.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 0.383u"mol/L"]]
)
=#
#=
reaction_order = 3
trials = Dict(
    250.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 1.454u"mol/L"]],
    300.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 1.258u"mol/L"]],
    350.0u"°C" => [[0.0u"s", 2.0u"mol/L"], [60.0u"s", 0.812u"mol/L"]]
)
=#
T_k_pairs = Dict()

for (temperature, time_and_concentrations) in trials
    time_vals = [time_and_concentration[1] for time_and_concentration in time_and_concentrations]
    concentration_vals = [time_and_concentration[2] for time_and_concentration in time_and_concentrations]

    if reaction_order == 0
        T_k_pairs[temperature] = -(concentration_vals[2] - concentration_vals[1]) / (time_vals[2] - time_vals[1])
    elseif reaction_order == 1
        T_k_pairs[temperature] = -(log(ustrip(concentration_vals[2])) - log(ustrip(concentration_vals[1]))) / (time_vals[2] - time_vals[1])
    else #reaction_order == 2
        T_k_pairs[temperature] = (1/(concentration_vals[2]) - 1/(concentration_vals[1])) / (time_vals[2] - time_vals[1])
    end
end

T_k_pairs

one_over_temperatures = Float64[]
ln_rate_constants = Float64[]
"""
for (temperature, rate_constant) in T_k_pairs
    push!(one_over_temperatures, 1.0/ustrip(uconvert(u"K", temperature)))
    push!(ln_rate_constants, log(ustrip(rate_constant)))
end
"""

one_over_temperatures
ln_rate_constants



plt = plot(one_over_temperatures, ln_rate_constants,
    seriestype=:scatter,
    xlabel="1/T (K⁻¹)",
    ylabel="ln(k)",
    title="Arrhenius Plot for Esterification",
    legend=false,
    framestyle=:box,
    #xlims = (0, one_over_temperatures[end]) #if you want to see how you would get that y intercept
)

data = DataFrame(
    one_over_T = one_over_temperatures,
    ln_k = ln_rate_constants,
)

model = lm(@formula(ln_k ~ one_over_T), data)

coefficients = coef(model)

y_intercept = coefficients[1]  #this is ln(pre exponential factor)
pre_exponential_factor = exp(y_intercept)

slope = coefficients[2] #this is our Ea/R_GAS
Ea = slope / 8.314

