using Revise
using FVMFramework
using Pkg

using Ferrite
using DifferentialEquations
using SparseArrays
using ComponentArrays
using NonlinearSolve
import SparseConnectivityTracer, ADTypes
using ILUZero

grid_dimensions = (1, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

addcellset!(grid, "internal_cells", x -> x != "chicken") # all

# methanol, water, carbon_monoxide, hydrogen, carbon_dioxide
initial_mass_fractions = [1.0, 1.3, 0.0001, 0.0001, 0.0001]

initial_mass_fractions = initial_mass_fractions ./ sum(initial_mass_fractions)

van_t_hoff_A_vec = [5.39e-4, 1.50e-5, 1.05e-11]
van_t_hoff_dH_vec = [-46800.0, -115000.0, -110000.0]

MSR_rxn = MSRReaction(
    49500.0,    # Delta H (Endothermic) [J/mol]
    -3800.0,    # Delta Gibbs free at 298.15K [J/mol]
    298.15,     # ref temp [K]
    1.25e7,     # kf_A (Pre-exponential factor)
    103000.0,   # kf_Ea (Activation Energy) [J/mol]
    [1, 2],     # reactant_ids: Methanol, Water
    [1, 1],     # reactant_stoich_coeffs
    [5, 4],     # product_ids: CO2, Hydrogen
    [1, 3],     # product_stoich_coeffs: 1 CO2 + 3 H2
    [-1, -1, 0, 3, 1], # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)

MD_rxn = MDReaction(
    90200.0,    # Delta H (Endothermic) [J/mol]
    24800.0,    # Delta Gibbs free at 298.15K [J/mol]
    298.15,     # ref temp [K]
    1.15e11,    # kf_A
    170000.0,   # kf_Ea [J/mol]
    [1],        # reactant_ids: Methanol
    [1],        # reactant_stoich_coeffs
    [3, 4],     # product_ids: CO, Hydrogen
    [1, 2],     # product_stoich_coeffs: 1 CO + 2 H2
    [-1, 0, 1, 2, 0], # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)

WGS_rxn = WGSReaction(
    -41100.0,   # Delta H (Exothermic) [J/mol]
    -28600.0,   # Delta Gibbs free at 298.15K [J/mol]
    298.15,     # ref temp [K]
    3.65e7,     # kf_A (Note: often adjusted depending on specific catalyst)
    87500.0,    # kf_Ea [J/mol]
    [3, 2],     # reactant_ids: CO, Water
    [1, 1],     # reactant_stoich_coeffs
    [5, 4],     # product_ids: CO2, Hydrogen
    [1, 1],     # product_stoich_coeffs
    [0, -1, -1, 1, 1], # stoich coefficients: [MeOH, H2O, CO, H2, CO2]
    van_t_hoff_A_vec,  # A vector (CH3O, HCOO, OH)
    van_t_hoff_dH_vec # dH vector (CH3O, HCOO, OH) [J/mol]
)

reaction_physics = SimpleChemPhysics(0.8, 4184, [MSR_rxn, MD_rxn, WGS_rxn], [1250, 1250, 1250], [0], 0)

phys_vec = [reaction_physics]

internal_cell_set_idxs = Set(getcellset(grid, "internal_cells"))

struct CellSet
    props_id::Int
    cell_set_idxs::Set{Int}
end

internal_cell_set = CellSet(1, internal_cell_set_idxs)

cell_sets = [internal_cell_set]

free_idxs = Int[]

temp_fixed_idxs = Int[]
chem_fixed_idxs = Int[]
#should probably create separate vectors for situations where we need to fix temperature but not mass fractions

function my_bc_mapper(cell_id)
    if cell_id in cell_sets[1].cell_set_idxs
        chem_bc = ChemBC(initial_mass_fractions)
        heat_bc = HeatBC(553.0)
        push!(temp_fixed_idxs, cell_id)
        push!(free_idxs, cell_id)
        return [chem_bc, heat_bc]
    end
end

heat_bcs = HeatBC[]
chem_bcs = ChemBC[]

n_cells = length(grid.cells)

for cell_id in 1:n_cells
    bcs = my_bc_mapper(cell_id) #returns vector [BC1, BC2]
    push!(heat_bcs, bcs[2])
    push!(chem_bcs, bcs[1])
end

boundary_map = SimpleReactionPhysicsBCs(heat_bcs, chem_bcs)

bc_sys = SimpleReactionBoundarySystem(boundary_map, free_idxs, temp_fixed_idxs, chem_fixed_idxs)

n_reactions = length(reaction_physics.chemical_reactions)
n_species = length(initial_mass_fractions)
alloc_mass_fraction_vec = zeros(n_species, n_cells)

u_proto = ComponentArray(mass_fractions=alloc_mass_fraction_vec, temp=zeros(n_cells))

for cell_id in 1:n_cells
    u_proto.mass_fractions[:, cell_id] = chem_bcs[cell_id].initial_mass_fractions
    u_proto.temp[cell_id] = heat_bcs[cell_id].initial_temp
end

u_axes = getaxes(u_proto)[1]
#we use u_proto and u_axes in the function because KrylovJL_GMRES complains when a component array from ComponentArrays.jl is passed in

u0 = Vector(u_proto)

du0 = u0 .* 0.0

cell_props_id_map = Int[]

for cell in CellIterator(grid)
    cell_id = cellid(cell)
    for i in eachindex(cell_sets)
        if cell_id in cell_sets[i].cell_set_idxs
            push!(cell_props_id_map, cell_sets[i].props_id)
        end
    end
end

initial_node_coordinates = get_node_coordinates(grid)

cell_neighbor_map, neighbor_map_respective_node_ids = get_neighbor_map(grid)

unconnected_cell_face_map, unconnected_map_respective_node_ids = get_unconnected_map(grid)

nodes_of_cells = get_nodes_of_cells(grid)

cell_volumes,
cell_centroids, #cell volumes and cell centroids are accessed at the id of the cell
connection_areas,
connection_normals,
connection_distances, #connection areas, normals, and distances are simply accessed by their location in the list which corresponds to the respective connection in cell_neighbor_map
unconnected_areas,
unconnected_normals = rebuild_fvm_geometry(
    cell_neighbor_map, neighbor_map_respective_node_ids,
    unconnected_cell_face_map, unconnected_map_respective_node_ids,
    initial_node_coordinates, nodes_of_cells
)

molar_concentrations_cache = zeros(Float64, length(u_proto.mass_fractions[:, 1])) #just using mass fractions for cell 1, this may cause some issues later!
net_rates_cache = zeros(Float64, length(reaction_physics.chemical_reactions))
#passing in caches into the FVM_iter_f! function still seems to be impossible, I guess we'll need PreallocationTools

f_closure_implicit = (du, u, p, t) -> simple_reaction_0D_f!(
    du, u, p, t,
    cell_neighbor_map,
    cell_volumes, cell_centroids,
    connection_areas, connection_normals, connection_distances,
    unconnected_areas,
    [0.03204, 0.01802, 0.02801, 0.00202, 0.04401],
    cell_props_id_map, bc_sys, phys_vec,
    u_axes, n_reactions, n_species
)
#just re-add t to the FVM_iter_f! function above to make it compatible with implicit solving

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

p_guess = 0.0

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0, u0, detector)

ode_func = ODEFunction(f_closure_implicit, jac_prototype=float.(jac_sparsity))

t0, tMax = 0.0, 1000000000.0
desired_steps = 10
dt = tMax / desired_steps
tspan = (t0, tMax)
t = t0:dt:tMax;

implicit_prob = ODEProblem(ode_func, u0, tspan, p_guess)

@time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true))

methanol_mass_fractions = [u_vec[1] for u_vec in sol.u]
methanol_initial = methanol_mass_fractions[1]
methanol_final = methanol_mass_fractions[end]
methanol_conversion = 1 - methanol_final / methanol_initial

#@VSCodeServer.profview sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))

#=
explicit_prob = ODEProblem(f_closure_implicit, u0, tspan, p_guess)
sol = solve(explicit_prob, Tsit5())
=#