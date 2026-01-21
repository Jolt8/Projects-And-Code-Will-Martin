using Revise
using Pkg
#Pkg.develop(path="c:/Users/wille/Desktop/Legacy-Projects-And-Code-Will-Martin/FVMFramework") 

using FVMFramework

using Ferrite
using DifferentialEquations
using SparseArrays
using ComponentArrays
using NonlinearSolve
import SparseConnectivityTracer, ADTypes
using ILUZero

grid_dimensions = (2, 2, 2)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

addcellset!(grid, "internal_cells", x -> x != "chicken") # all

# acetic acid, ethanol, acetic_acid, water
initial_mass_fractions = [0.5, 0.5, 0.0, 0.0]

methanol_reforming_rxn = ChemicalReaction(
    -6.0, #Delta H, -1.0 = exothermic
    -4.0, #Delta Gibbs free at ref temp
    298, #ref temp
    1000, #kf_A
    62, #kf_Ea
    [1, 2], #reactant_ids
    [1, 1], #reactant_stoich_coeffs
    [3, 4], #product_ids
    [1, 1], #product_stoich_coeffs
    [-1, -1, 1, 1] #stoich coefficients -1 = reactant, 1 = product
)
#=
rwgs_rxn = ChemicalReaction(
    #stuff here later
)
=#

reaction_physics = ChemPhysics(0.6e-3, 1000, 4.184, [methanol_reforming_rxn], [0], 0) 

chem_phys_vec = [reaction_physics]

internal_cell_set_idxs= Set(getcellset(grid, "internal_cells"))

struct CellSet
    props_id::Int
    cell_set_idxs::Set{Int}
end

internal_cell_set = CellSet(1, internal_cell_set_idxs)

cell_sets = [internal_cell_set]

free_idxs = Int[]
dirichlet_idxs = Int[]
#should probably create separate vectors for situations where we need to fix temperature but not mass fractions

function my_bc_mapper(cell_id)
    if cell_id in cell_sets[1].cell_set_idxs
        chem_bc = ChemBC(initial_mass_fractions)
        heat_bc = HeatBC(500.0)
        push!(free_idxs, cell_id)
        return [chem_bc, heat_bc] 
    end
end

chem_bcs = ChemBC[]
heat_bcs = HeatBC[]

n_cells = length(grid.cells)

for cell_id in 1:n_cells
    bcs = my_bc_mapper(cell_id) #returns vector [BC1, BC2]
    push!(chem_bcs, bcs[1])
    push!(heat_bcs, bcs[2])
end

boundary_map = MultiPhysicsBCs(chem_bcs, heat_bcs)

bc_sys = BoundarySystem(boundary_map, free_idxs, dirichlet_idxs)

n_reactions = length(reaction_physics.chemical_reactions)
n_species = length(initial_mass_fractions)
alloc_mass_fraction_vec = zeros(n_species, n_cells) 

u_proto = ComponentArray(mass_fractions = alloc_mass_fraction_vec, temp = zeros(n_cells))

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

f_closure_implicit = (du, u, p, t) -> methanol_reformer_f!(
    du, u, p, t,
    cell_neighbor_map,
    cell_volumes, cell_centroids, 
    connection_areas, connection_normals, connection_distances, 
    unconnected_areas, 
    [0.06005, 0.04607, 0.08811, 0.018015],
    cell_props_id_map, bc_sys, chem_phys_vec, 
    molar_concentrations_cache, net_rates_cache,
    u_axes, n_reactions, n_species
)
#just re-add t to the FVM_iter_f! function above to make it compatible with implicit solving

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

p_guess = 0.0

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0, u0, detector)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0, tMax = 0.0, 10000000000.0 
desired_steps = 10
dt = tMax / desired_steps
tspan = (t0, tMax)
t = t0:dt:tMax;

implicit_prob = ODEProblem(ode_func, u0, tspan, p_guess)

@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
@VSCodeServer.profview sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))