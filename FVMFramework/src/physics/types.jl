### Physics ###
abstract type AbstractPhysics end

## Chemistry Physics ##
struct ChemicalReaction
    heat_of_reaction::Float64
    delta_gibbs_free_energy::Float64
    K_gibbs_free_ref_temp::Float64
    kf_A::Float64 # Pre-exponential factor for forward reaction
    kf_Ea::Float64 # Activation energy for forward reaction
    reactants::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    reactant_stoich_coeffs::Vector{Int} #(stoich_coeff_for_reactant_1, stoich_coeff_for_reactant_2, etc...)
    products::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    product_stoich_coeffs::Vector{Int} #[stoich_coeff_for_product_1, stoich_coeff_for_product_2, etc...]
    all_stoich_coeffs::Vector{Int} #(stoich_for_chemical_id_1, stoich_for_chemical_id_2, etc...) #-1 = reactant, 1 = product
end

struct ChemPhysics <: AbstractPhysics
    k::Float64
    rho::Float64
    cp::Float64
    chemical_reactions::Vector{ChemicalReaction}
    chemical_vol_source_term::Vector{Float64} #chemical addition
    heat_vol_source_term::Float64 #volumetric heating
end

## Heat Physics ##
struct HeatPhysics <: AbstractPhysics
    k::Float64 
    rho::Float64
    cp::Float64
    heat_vol_source_term::Float64 #volumetric heating
end

struct FluidPhysics <: AbstractPhysics
    k::Float64 
    rho::Float64
    mu::Float64
    cp::Float64
end




### BC TYPES ###

abstract type AbstractBC end

## Chemistry BCs ##
struct ChemBC <: AbstractBC
    initial_mass_fractions::Vector{Float64}
end

## Heat BCs ##
struct HeatBC <: AbstractBC
    initial_temp::Float64
end

struct VelBC <: AbstractBC
    type::Symbol 
    initial::Float64
end

struct PressureBC <: AbstractBC
    type::Symbol 
    initial::Float64
end

## General BC Types ## 

struct MultiPhysicsBCs #this also defines the order of each variable in u used in the future
    vel_x_bcs::Vector{VelBC}
    vel_y_bcs::Vector{VelBC}
    vel_z_bcs::Vector{VelBC}
    pressure_bcs::Vector{PressureBC}
    temp_bcs::Vector{HeatBC}
    chem_bcs::Vector{ChemBC}
end

struct BoundarySystem
    boundary_map::MultiPhysicsBCs
    free_idxs::Vector{Int}
    dirichlet_idxs::Vector{Int}
end
