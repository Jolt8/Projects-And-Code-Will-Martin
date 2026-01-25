module FVMFramework

using Ferrite
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using SciMLSensitivity
using Optimization
using OptimizationPolyalgorithms
using Zygote
#using Enzyme
using RecursiveArrayTools
using OptimizationOptimJL
using ILUZero
using NonlinearSolve
using ComponentArrays
using StaticArrays
using ProfileView

import AlgebraicMultigrid
import SparseConnectivityTracer
import ADTypes 
import Logging

#this is here just in case we have to do this again
#add Ferrite DifferentialEquations LinearAlgebra SparseArrays SciMLSensitivity Optimization OptimizationPolyalgorithms Zygote Enzyme RecursiveArrayTools OptimizationOptimJL ILUZero NonlinearSolve ComponentArrays StaticArrays ProfileView
#add AlgebraicMultigrid SparseConnectivityTracer ADTypes Logging

#to add new files:
#create a terminal on this file
#type julia
#using Pkg
#Pkg.activate(".")
#then ] and add whatever package you'd like

include("geometry.jl")
export get_node_coordinates, get_neighbor_map, get_unconnected_map
export rebuild_fvm_geometry, get_nodes_of_cells

# Physics
include("physics/types.jl")
export AbstractPhysics, AbstractReaction, AbstractBoundarySystem
export SimpleChemPhysics, ChemPhysics
export HeatPhysics
export ChemBC, HeatBC, MultiPhysicsBCs

include("physics/helper_functions.jl")
export R_gas, upwind, harmonic_mean, van_t_hoff, arrenhius_k
export get_mw_avg, cell_rho_ideal, get_cell_cp

include("physics/advection.jl")
export species_advection!, enthalpy_advection!

include("physics/diffusion.jl")
export species_numerical_flux, diffusion_mass_fraction_exchange!

include("physics/darcy_flow.jl")
export get_darcy_mass_flux, continuity_and_momentum_darcy

include("physics/heat_transfer.jl")
export get_k_effective, numerical_flux, diffusion_temp_exchange!

include("physics/chemistry.jl")
export PowerLawReaction
export net_reaction_rate, K_gibbs_free, react_cell!

include("physics/methanol_reforming_net_rates.jl")
export MSRReaction, MDReaction, WGSReaction #new reaction types
export net_reaction_rate

include("solvers/preconditioners.jl")
export iluzero, algebraicmultigrid

include("solvers/fvm_operators/methanol_reformer_operator.jl")
export methanol_reformer_f!

include("solvers/fvm_operators/simple_reaction_0D.jl")
export simple_reaction_0D_f!, SimpleReactionBoundarySystem, SimpleReactionPhysicsBCs

end