### Physics ###
abstract type AbstractPhysics end

## Chemistry Physics ##
abstract type AbstractReaction end

struct SimpleChemPhysics <: AbstractPhysics
    rho::Float64
    cp::Float64
    chemical_reactions::Vector{AbstractReaction}
    cell_kg_cat_per_m3_for_each_reaction::Vector{Float64}
    chemical_vol_source_term::Vector{Float64} #chemical addition
    heat_vol_source_term::Float64 #volumetric heating
end

struct ChemPhysics <: AbstractPhysics
    k::Float64
    rho::Float64
    cp::Float64
    permeability::Float64
    chemical_reactions::Vector{AbstractReaction}
    cell_kg_cat_per_m3_for_each_reaction::Vector{Float64}
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

abstract type  MultiPhysicsBCs end 

abstract type AbstractBoundarySystem end
