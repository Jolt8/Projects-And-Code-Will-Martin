using ModelingToolkit
using DifferentialEquations
#using Plots
#using GLMakie


import ModelingToolkit: t_nounits as t, D_nounits as D
"""
using Clapeyron
model = PR(["water"])
test = @register_symbolic Clapeyron.volume(p, T, z, phase)
"""


#Maybe an effective method of generating the K_eq expression would be as follows 
function K_eq(forward_k, reverse_k, reactants_partial_pressures_vec, reactants_stoichs_vec, products_molar_flows_vec, products_stoichs_vec,)
    #reactants_stoichs_vec would be formatted [1]
    #reactants_partial_pressures_vec would be formatted [5900]
    #products_partial_pressures_vec would be formatted [11009, 32098]
    #products_stoichs_vec would be formatted [1, 2]
    #The only reason we're doing this is because dicts cannot be handled within mtkconnectors
    forward_term = 1
    for i in eachindex(reactants_stoichs_vec)
        forward_term *= (reactants_partial_pressures_vec[i]^reactants_stoichs_vec[i])
    end
    reverse_term = 1
    for i in eachindex(reactants_stoichs_vec)
        reverse_term *= (products_partial_pressures_vec[i]^products_stoichs_vec[i])
    end
    overall_rate = forward_k * forward_term - reverse_k * reverse_term
    return overall_rate
end
@register_symbolic K_eq()


species_val = ["CH3OH", "H2O", "CO", "H2", "CO2"]
n_species_val = length(species_val)

heat_rxn_vec_val = [-101, 54]


#Start of Component Definitions
@connector ReactorPort begin
    @parameters begin
        species[1:n_species_val] = [], [description = "species present in mixture"]
        F_guess[1:n_species_val] = [], [description = "molar flows in mol/s for each respective species"]
        T_guess = 273.13, [description = "temperature in K"]
        p_guess = 101325, [description = "pressure in Pa"]
    end

    @variables begin
        F(t)[1:n_species_val], [connect = Flow, guess = F_guess]
        T(t), [guess = T_guess]
        p(t), [guess = p_guess]
    end
end

@mtkmodel ReactionKinetics begin
    @parameters begin
        Ea_f
        Ea_r
        kf_A
        kf_Ea
        kr_A
        kr_Ea
        reactants_stoichs_vec, [description=""]
        products_stoichs_vec, [description=""]
        R_gas = 8.314
    end
    
    @variables begin
        F(t)
        T(t)
        p(t)
        rate(t)
        reactants_partial_pressures(t)[1:length(reactants_stoichs_vec)]
        products_partial_pressures(t)[1:length(products_stoichs_vec)]
        forward_terms(t)[1:length(reactants_stoichs_vec)]
        reverse_terms(t)[1:length(products_stoichs_vec)]
    end
    
    @equations  begin
        kf = kf_A * exp(-kf_Ea / (R_gas * T))
        kr = kr_A * exp(-kr_Ea / (R_gas * T))
        forward_terms[1:length(reactants_stoichs_vec)] *= [reactants_partial_pressures[i]^reactants_stoichs_vec[i] for i in eachindex(reactants_stoichs_vec)]
        reverse_terms[1:length(products_stoichs_vec)] *= [products_partial_pressures[i]^products_stoichs_vec[i] for i in eachindex(product_stoichs_vec)]
        rate ~ kf .* forward_term - kr .* reverse_term
    end
end

@mtkmodel ErgunPressureDrop begin
    @parameters begin
        G = 0.002, [description="Factor to relate total flow to mass velocity"]
        ϵ = 0.3,     [description="Bed void fraction"]
        μ = 2.0e-5,   [description="Viscosity Pa*s"]
        Dp = 0.005,    [description="Particle diameter m"]
        ρ = 1000.0,[description="Catalyst density kg/m^3"]
    end
    
    @variables begin
        p(t)
    end

    @equations begin
        term_1 = (150 * μ * G*(1 - ϵ)^2) / (ρ * Dp^2 * ϵ^3)
        term_2 = (1.75 * G^2 * (1 - ϵ)) / (ρ * Dp * ϵ^3)
        D(p) = -1 * (term_1 + term_2)
    end
end

#Note sure if a reactor segment requires a 1D Element to connect them 
#Perhaps this would be a good use case for trying out the finite volume method in mtkcomponents
@mtkmodel ReactorConnector begin
    @components begin
        inlet = ReactorPort()
        oulet = ReactorPort()
    end
end

@mtkmodel ReactorSegment begin
    @components begin
        port = ReactorPort()
        PressureDrop = ErgunPressureDrop()
    end

    @parameters begin
        catalyst_weight = 1.0, [description = "weight of catalyst in kg present in reactor segment"] 
        #I don't know if there's a method to convert catalyst weight in each segment to reactor length like ρ_bulk * cross_sect_area
        heat_of_reaction = -92.3
        heat_capacity_vec[1:n_species_val], [description = "[J/(g*K)]"]
        reactants_stoichs_vec[1:2]
        products_stoichs_vec[1:3]
    end

    @variables begin
        F(t)
        T(t) = 293.13
        p(t) = 101325
        der_T(t) = 0.0
        der_F(t)[1:n_species_val]
    end

    @equations begin
        T ~ port.T
        p ~ port.p
        D(F)[1:length(reactants_stoichs_vec)] ~ [-1 * rate * F(t)[i] for i in 1:length(reactants_stoichs_vec)] #Don't know how do a negative rate for reactants
        D(F)[length(reactants_stoichs_vec):length(products_stoichs_vec)] ~ [rate * F(t)[i] for i in length(reactants_stoichs_vec):length(products_stoichs_vec)]
        D(T) ~ (rate * -heat_of_reaction) / prod([F(t)[i] * heat_capacity_vec[i] for i in eachindex(F(t))])
        D(P) ~ PressureDrop.D(P)
    end
end

@mtkmodel SpeciesSource begin
    @components begin
        port = ReactorPort()
    end

    @parameters begin
        F_fixed[1:n_species] = [], [description = "molar flows in mol/s for each respective species"]
        T_fixed = 273.13, [description = "temperature in K"]
        p_fixed = 101325, [description = "pressure in Pa"]
    end

    @equations begin
        port.F ~ F_fixed
        port.T ~ T_fixed
        port.p ~ p_fixed
    end
end

rows = 4
n_nodes = rows

heat_capacity_vec_val = [20, 20, 20, 20, 20]
reactants_stoichs_vec_val = [1]
products_stoichs_vec_val = [1, 2]

nodes = [ReactorSegment(name=Symbol("Reactor_Segment_", i), catalyst_weight=1, heat_of_reaction=-92.0, heat_capacity_vec=heat_capacity_vec_val, reactants_stoichs_vec=reactants_stoichs_vec_val, products_stoichs_vec=products_stoichs_vec_val) for i in 1:rows]

# Add boundary condition components
left_bcs = SpeciesSource(F_fixed=[1.0, 1.3, 0.001, 0.001, 0.001], T_fixed=293.15, p_fixed=101325, name=Symbol("Species_Source"))
right_bcs = SpeciesSource(F_fixed=[1.0, 1.3, 0.001, 0.001, 0.001], T_fixed=293.15, p_fixed=101325, name=Symbol("Species_Source"))

connections = Equation[]

push!(connections, connect(left_bcs.port, nodes[1].port))
[push!(connections, connect(nodes[i].port, nodes[i].port)) for i in 1:(rows-1)]
push!(connections, connect(nodes[end].port, right_bcs.port))

# --- Simulation Setup ---
eqs_total = System[]
all_systems = vcat(
    vec(nodes),
    vec(left_bcs),
    vec(right_bcs)
)
all_systems

@named rod = ODESystem(connections, t, systems=all_systems)

sys = structural_simplify(rod)  

tspan = (0.0, 10.0)

prob = ODEProblem(sys, [], tspan)

sol = solve(prob)
