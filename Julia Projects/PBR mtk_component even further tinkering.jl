using ModelingToolkit
using DifferentialEquations
#using Plots
#using GLMakie


import ModelingToolkit: t_nounits as t, D_nounits as D
#Start of Component Definitions
@connector ReactorPort begin
    @structural_parameters begin
        n_species
        n_reactants
        n_products
    end

    @parameters begin
        F_guess[1:n_species], [connect = Flow]
        T_guess = 273.13, [description = "temperature in K"]
        p_guess = 101325, [description = "pressure in Pa"]
    end

    @variables begin
        F(t)[1:n_species] = [1.0 for i in 1:n_species] #, [connect = Flow, guess = F_guess]
        T(t), [guess = T_guess]
        p(t), [guess = p_guess]
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
        D(p) ~ -1 * ((150 * μ * G*(1 - ϵ)^2) / (Dp^2 * ϵ^3) + (1.75 * G^2 * (1 - ϵ)) / (Dp * ϵ^3))
    end
end

#Note sure if a reactor segment requires a 1D Element to connect them 
#Perhaps this would be a good use case for trying out the finite volume method in mtkcomponents
@mtkmodel ReactorSegment begin
    @structural_parameters begin
        n_species
        n_reactants
        n_products
    end

    @variables begin
        F(t)[1:n_species] = [1.0 for i in 1:n_species]
        T(t) = 293.13
        p(t) = 101325

        #kinetics part
        F_total(t) = 1.0

        kf(t) = 1.0
        kr(t) = 1.0

        partial_pressures(t)[1:n_species] = [1.0 for i in 1:n_species]

        pre_f_t(t)[1:n_species] = [1.0 for i in 1:n_species]
        pre_r_t(t)[1:n_species] = [1.0 for i in 1:n_species]
        forward_term(t) = 1.0
        reverse_term(t) = 1.0

        rate(t) = 1.0
    end

    @parameters begin
        catalyst_weight = 1.0, [description = "weight of catalyst in kg present in reactor segment"] 
        heat_of_reaction = -92.3
        heat_capacity_vec[1:n_species] = [1.0 for i in 1:n_species]
        stoich_coeffs[1:n_species] = [1.0 for i in 1:n_species]
        
        reactants_stoichs[1:n_species] = [1.0 for i in 1:n_species]
        products_stoichs[1:n_species] = [1.0 for i in 1:n_species]
        kf_A = 1.0
        kf_Ea = 1.0
        kr_A = 1.0
        kr_Ea = 1.0

        R_gas = 8.134

        F_in[1:n_species] = [1.0 for i in 1:n_species]
        T_in = 293.13
        p_in = 101325.0
    end

    @equations begin
        F_total ~ sum(F)
        
        kf ~ kf_A * exp(-kf_Ea / (R_gas * T))
        kr ~ kr_A * exp(-kr_Ea / (R_gas * T))
        
        partial_pressures[1:n_species] ~ [p * F[i] / F_total for i in 1:n_species] #this is the only syntax that works for doing arrays within parameters
        
        pre_f_t[1:n_species] ~ [partial_pressures[i]^reactants_stoichs[i] for i in 1:n_species] #this is a pretty brittle way because we're relying on num^0 = 1, Too Bad!
        pre_r_t[1:n_species] ~ [partial_pressures[i]^products_stoichs[i] for i in 1:n_species]
        forward_term ~ 10 #prod(pre_f_t)
        reverse_term ~ 10 #prod(pre_r_t)
        #NOTE: WHAT THE FUCK: Why does the forward and reverse only work if we do a "pre" version of them
        #perhaps reduce(*, pre_f_t)

        rate ~ kf * forward_term - kr * reverse_term
        #pressure_drop.p ~ p

        D(F) ~ [F_in[i] - F[i] + rate * catalyst_weight * stoich_coeffs[i] for i in 1:n_species] #don't use dot here: [rate .* stoich_coeffs]

        D(T) ~ ((sum(F_in .* heat_capacity_vec) * T_in - sum(F .* heat_capacity_vec) * T) + 
                rate * catalyst_weight * (-heat_of_reaction)) / 
                sum(F .* heat_capacity_vec) #wait what, why does dot work here?
        #NOTE: getting a cannot 'convert' error message usually means you forgot to use ~
    end
end


@mtkmodel SpeciesSource begin
    @components begin
        port = ReactorPort(n_species=n_species, n_reactants=n_reactants, n_products=n_products)
    end

    @structural_parameters begin
        n_species
        n_reactants
        n_products
    end

    @parameters begin
        F_fixed[1:n_species] = [1.0 for i in 1:n_species], [connect = Flow]
        T_fixed = 273.13, [description = "temperature in K"]
        p_fixed = 101325, [description = "pressure in Pa"]
    end

    @equations begin 
        port.F ~ F_fixed
        port.T ~ T_fixed
        port.p ~ p_fixed
    end
end

stoich_coeffs_val = [-1, 0, 1, 2, 0] # Negative for reactants, positive for products

n_species_val = length(stoich_coeffs_val)
n_reactants_val = 1
n_products_val = 2

rows = 1
n_nodes = rows

reactants_stoichs_val = [1, 0, 0, 0, 0]
products_stoichs_val = [0, 0, 1, 2, 0]

@named left_bcs = SpeciesSource(n_species=n_species_val, n_reactants=n_reactants_val, n_products=n_products_val, F_fixed=[1.0, 1.0, 0.01, 0.01, 0.01], T_fixed=500.0, p_fixed=5e5)

#NOTE: whenever it says that it cannot convert a type into an equation, it usually means you are incorrectly defining equation as not an equation somewhere in your model
@named test = ReactorSegment(name=Symbol("seg"),
    n_species = n_species_val,
    n_reactants = 1,
    n_products = 1,
    catalyst_weight = 5,
    heat_of_reaction = 90.7e3, # J/mol
    heat_capacity_vec = [81.6, 33.6, 29.1, 28.8, 37.1], # J/(mol*K)
    stoich_coeffs = stoich_coeffs_val,
    reactants_stoichs=reactants_stoichs_val,
    products_stoichs=products_stoichs_val,
    kf_A=1.0e5, kf_Ea=8.0e4, kr_A=1.0e12, kr_Ea=1.5e5)

println(ModelingToolkit.get_defaults(test))
structural_simplify(test)


segments = [ReactorSegment(name=Symbol("seg", i),
    n_species = n_species_val,
    n_reactants = 1,
    n_products = 1,
    catalyst_weight = 5,
    heat_of_reaction = 90.7e3, # J/mol
    heat_capacity_vec = [81.6, 33.6, 29.1, 28.8, 37.1], # J/(mol*K)
    stoich_coeffs = stoich_coeffs_val,
    reactants_stoichs=reactants_stoichs_val,
    products_stoichs=products_stoichs_val,
    kf_A=1.0e5, kf_Ea=8.0e4, kr_A=1.0e12, kr_Ea=1.5e5) for i in 1:rows
]
#println(defaults(segments[1]))

#NOTE: in the future, if you want to change a value like kf_A halfway through the reactor just change kinetics directly like: seg.kinetics.kf_A

#nodes = [ReactorSegment(name=Symbol("Reactor_Segment_", i), catalyst_weight=1, heat_of_reaction=-92.0, heat_capacity_vec=heat_capacity_vec_val, reactants_stoichs=reactants_stoichs_val, products_stoichs_vec=products_stoichs_vec_val) for i in 1:rows]

#

#!!!!IMPORTANT!!!! NOTE: For some reason, whenever you do F(t)[1:n] YOU CANNOT DO A connect = Flow !!!!IMPORTANT!!!!
#actually, I'm not sure if the above is true, this requires more testing in a later complete system

@named right_bcs = SpeciesSource(n_species=n_species_val, n_reactants=n_reactants_val, n_products=n_products_val, F_fixed=[1.0, 1.0, 0.01, 0.01, 0.01], T_fixed=500.0, p_fixed=5e5)

connections = Equation[]
# Left boundary condition
for i in 1:n_species_val
    push!(connections, segments[1].F_in[i] ~ left_bcs.port.F[i])
end
push!(connections, segments[1].T_in ~ left_bcs.port.T)
push!(connections, segments[1].p_in ~ left_bcs.port.p)

# Inter-segment connections (outlet of i feeds inlet of i+1)
for i in 1:(rows-1)
    for j in 1:n_species_val
        push!(connections, segments[i+1].F_in[j] ~ segments[i].F[j])
    end
    #push!(connections, segments[i+1].F_in ~ segments[i].F)
    push!(connections, segments[i+1].T_in ~ segments[i].T)
    push!(connections, segments[i+1].p_in ~ segments[i].p)
end
"""
# Right boundary (just observation, could be removed)
push!(connections, right_bcs.port.F ~ segments[end].F)
push!(connections, right_bcs.port.T ~ segments[end].T)
push!(connections, right_bcs.port.p ~ segments[end].p)"""
# --- Simulation Setup ---
eqs_total = System[]
all_systems = vcat(
    segments,
    left_bcs
)
all_systems

@named rod = ODESystem(connections, t, systems=all_systems)

println(defaults(rod))

sys = structural_simplify(rod)

println(defaults(sys))

tspan = (0.0, 10.0)

prob = ODEProblem(sys, [], tspan)

sol = solve(prob)
