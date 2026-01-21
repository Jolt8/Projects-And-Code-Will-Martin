using Ferrite
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using SciMLSensitivity
using Optimization, OptimizationPolyalgorithms, Zygote
using Enzyme
using RecursiveArrayTools
using OptimizationOptimJL
using ILUZero
import AlgebraicMultigrid
import SparseConnectivityTracer, ADTypes 
using NonlinearSolve
import Logging
using ComponentArrays
using StaticArrays
using ProfileView

#add Ferrite DifferentialEquations LinearAlgebra SparseArrays SciMLSensitivity Optimization OptimizationPolyalgorithms Zygote Enzyme RecursiveArrayTools OptimizationOptimJL ILUZero NonlinearSolve ComponentArrays StaticArrays ProfileView

function get_node_coordinates(grid)
    n_nodes = length(grid.nodes)

    node_coordinates = Vector{SVector{3, Float64}}(undef, n_nodes)

    visited_map = Set{Int}()

    for node_id in 1:length(grid.nodes)
        if !(node_id in visited_map)
            #add it to the set so its coordinates do not get re-added
            push!(visited_map, node_id)
            coordinates = get_node_coordinate(grid.nodes[node_id])
            
            node_coordinates[node_id] = SVector(coordinates[1], coordinates[2], coordinates[3])
        end
    end
    return node_coordinates
end

function get_cell_topology(grid)
    n_nodes = length(grid.nodes)

    topology = Vector{Vector{Int}}(undef, n_nodes)

    visited_map = Set()

    for cell in CellIterator(grid)
        cell_id = cellid(cell)
        for i in 1:length(cell.nodes)
            node_id = cell.nodes[i]
            if !(node_id in visited_map)
                push!(visited_map, node_id) #add it to the set so its coordinates do not get re-added

                topology[node_id] = [cell_id]
            else
                push!(topology[node_id], cell_id)
            end
        end
    end
    return topology
end

function get_nodes_of_cells(grid)
    n_cells = length(grid.cells)

    node_topology = Vector{Vector{Int}}(undef, n_cells)

    for cell in CellIterator(grid)
        cell_id = cellid(cell)
        node_topology[cell_id] = collect(grid.cells[cell_id].nodes)
    end
    
    return node_topology
end

function get_face_nodes(grid, cell_id, face_idx)
    cell = grid.cells[cell_id]
    face_nodes = collect(Ferrite.faces(cell)[face_idx])
    return face_nodes
end

function get_cell_faces_node_ids(grid)
    n_cells = length(grid.cells)

    cell_faces_node_ids = NTuple{6, NTuple{4, Int64}}[]
    #Vector{NTuple{6, NTuple{4, Int64}}}

    for cell_id in 1:n_cells
        current_cell_faces_node_ids = Ferrite.faces(grid.cells[cell_id])
        push!(cell_faces_node_ids, current_cell_faces_node_ids)
    end
    return cell_faces_node_ids
end

#get_cell_faces_node_ids(grid)

function get_neighbor_map(grid)
    nodes_of_cells = get_nodes_of_cells(grid)
    #returns:
    #cell/neighbor pairs (ex. (1, 2) (cell_1 is connected with cell_2))
    #respective node idxs of cell and neighbor (ex. (2, 8, 44, 38))

    cell_neighbor_pairs = Tuple{Int, Int}[]
    neighbor_map_respective_node_ids = NTuple{4, Int}[]

    n_cells = length(grid.cells)

    top = ExclusiveTopology(grid)
    for cell_id in 1:n_cells
        for face_idx in 1:nfacets(grid.cells[cell_id])
            neighbor_info = top.face_face_neighbor[cell_id, face_idx]

            if !isempty(neighbor_info)
                neighbor_id = collect(neighbor_info[1].idx)[1]
                
                if cell_id < neighbor_id

                    curr_cell_nodes = nodes_of_cells[cell_id]
                    neighbor_nodes = nodes_of_cells[neighbor_id]

                    respective_nodes = get_face_nodes(grid, cell_id, face_idx)
                    
                    push!(cell_neighbor_pairs, (cell_id, neighbor_id))
                    push!(neighbor_map_respective_node_ids, ntuple(i -> respective_nodes[i], 4))
                end
            end
        end
    end
    return cell_neighbor_pairs, neighbor_map_respective_node_ids
end

function get_unconnected_map(grid)
    nodes_of_cells = get_nodes_of_cells(grid)
    #returns:
    #list of unconnected faces for each cell (ex. (1, 5) (cell_1's 5th face_idxs is not connected))
    #respective node idxs of cell face (ex. (2, 8, 44, 38))

    n_cells = length(grid.cells)

    unconnected_cell_face_map = NTuple{2, Int}[]
    unconnected_map_respective_node_ids = NTuple{4, Int}[]

    n_cells = length(grid.cells)

    top = ExclusiveTopology(grid)
    for cell_id in 1:n_cells
        for face_idx in 1:nfacets(grid.cells[cell_id])
            neighbor_info = top.face_face_neighbor[cell_id, face_idx]

            if isempty(neighbor_info) #note that were checking if it isempty here, not !isempty like above
                curr_cell_nodes = nodes_of_cells[cell_id]
                respective_nodes = get_face_nodes(grid, cell_id, face_idx)

                push!(unconnected_cell_face_map, (cell_id, face_idx))
                push!(unconnected_map_respective_node_ids, ntuple(i -> respective_nodes[i], 4))
            end
        end
    end
    return unconnected_cell_face_map, unconnected_map_respective_node_ids
end

function cross_product(a, b)
    x = a[2]*b[3] - a[3]*b[2]
    y = a[3]*b[1] - a[1]*b[3]
    z = a[1]*b[2] - a[2]*b[1]
    return SVector(x, y, z)
end

function calculate_hex_volume(p)
    c = sum(p) / 8.0
    
    # Faces of Hex8 (node indices for each face)
    faces = (
        (1, 4, 3, 2),  # bottom
        (1, 2, 6, 5),  # front
        (2, 3, 7, 6),  # right
        (3, 4, 8, 7),  # back
        (4, 1, 5, 8),  # left
        (5, 6, 7, 8)   # top
    )
    
    total_vol = 0.0
    for face in faces
        node_1, node_2, node_3, node_4 = p[face[1]], p[face[2]], p[face[3]], p[face[4]]
        
        # Split quad face into 2 triangles, form tetrahedra with centroid
        total_vol += dot(node_1 - c, cross_product(node_2 - c, node_3 - c))
        total_vol += dot(node_1 - c, cross_product(node_3 - c, node_4 - c))
    end
    
    return abs(total_vol) / 6.0
end

function calculate_face_area(face_node_coordiantes)
    node_1_coords = face_node_coordiantes[1]
    node_2_coords = face_node_coordiantes[2]
    node_3_coords = face_node_coordiantes[3]
    node_4_coords = face_node_coordiantes[4]

    cross_a = cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)
    cross_b = cross_product(node_3_coords - node_1_coords, node_4_coords - node_1_coords)
    
    area_vec_1 = 0.5 * cross_a
    area_vec_2 = 0.5 * cross_b

    total_area_vec = area_vec_1 + area_vec_2
    return norm(total_area_vec)
end

function rebuild_fvm_geometry(
        cell_neighbor_map, neighbor_map_respective_node_ids, 
        unconnected_cell_face_map, unconnected_map_respective_node_ids, 
        node_coordinates, nodes_of_cells
    )
    #The cells_data and connections map are still causing GC
    CoordType = eltype(node_coordinates)
    T = eltype(CoordType)

    n_cells = length(nodes_of_cells)

    cell_volumes = Vector{T}(undef, n_cells)
    cell_centroids = Vector{CoordType}(undef, n_cells)

    for cell_id in eachindex(nodes_of_cells)
        cell_nodes = nodes_of_cells[cell_id]

        n_nodes = length(cell_nodes)

        p = ntuple(8) do i
            @inbounds node_coordinates[cell_nodes[i]]
        end

        vol = calculate_hex_volume(p)
        
        cent = sum(p) / T(n_nodes)

        cell_volumes[cell_id] = vol
        cell_centroids[cell_id] = cent
    end

    n_connections = length(cell_neighbor_map)
    
    connection_areas = Vector{T}(undef, n_connections)
    connection_normals = Vector{CoordType}(undef, n_connections)
    #the above connection_normals causes some GC (4% of optimization runtime)
    connection_distances = Vector{T}(undef, n_connections)
    

    for (i, (cell_id, neighbor_id)) in enumerate(cell_neighbor_map)
        face_node_indices = neighbor_map_respective_node_ids[i] #neighbor_map_respective_node_ids[i] looks like (1, 4, 7, 21) 
        node_1_coords = node_coordinates[face_node_indices[1]]
        node_2_coords = node_coordinates[face_node_indices[2]]
        node_3_coords = node_coordinates[face_node_indices[3]]
        node_4_coords = node_coordinates[face_node_indices[4]]

        #get_area
        cross_a = cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)
        cross_b = cross_product(node_3_coords - node_1_coords, node_4_coords - node_1_coords)
        
        area_vec_1 = 0.5 * cross_a
        area_vec_2 = 0.5 * cross_b

        total_area_vec = area_vec_1 + area_vec_2
        total_area = norm(total_area_vec)

        #get normal
        cell_normal = normalize(total_area_vec)

        #get distance 
        dist = norm(cell_centroids[cell_id] - cell_centroids[neighbor_id])
        
        connection_areas[i] = total_area
        connection_normals[i] = cell_normal
        connection_distances[i] = dist
    end

    n_unconnected_faces = length(unconnected_cell_face_map)
    
    unconnected_areas = Vector{SVector{6, T}}(undef, n_unconnected_faces)
    unconnected_normals = Vector{SVector{6, CoordType}}(undef, n_unconnected_faces)

    for (i, (cell_id, face_idx)) in enumerate(unconnected_cell_face_map)
        areas = MVector{6, T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        normals = MVector{6, CoordType}(
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0), 
            (0.0, 0.0, 0.0)
        )
        for face_idx in 1:6
            face_node_indices = unconnected_map_respective_node_ids[face_idx] #unconnected_map_respective_node_ids[i] looks like (1, 4, 7, 21) 
            node_1_coords = node_coordinates[face_node_indices[1]]
            node_2_coords = node_coordinates[face_node_indices[2]]
            node_3_coords = node_coordinates[face_node_indices[3]]
            node_4_coords = node_coordinates[face_node_indices[4]]

            #get_area
            cross_a = cross_product(node_2_coords - node_1_coords, node_3_coords - node_1_coords)
            cross_b = cross_product(node_3_coords - node_1_coords, node_4_coords - node_1_coords)
            
            area_vec_1 = 0.5 * cross_a
            area_vec_2 = 0.5 * cross_b

            total_area_vec = area_vec_1 + area_vec_2
            total_area = norm(total_area_vec)

            areas[face_idx] = total_area 

            normals[face_idx] = normalize(total_area_vec)
        end
        unconnected_areas[cell_id] = areas

        unconnected_normals[cell_id] = normals
    end
    return cell_volumes, cell_centroids, connection_areas, connection_normals, connection_distances, unconnected_areas, unconnected_normals
end

abstract type AbstractPhysics end

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
    all_stoich_coeffs::Vector{Int} #(stoich_for_chemical_id_1, stoich_for_chemical_id_2, etc...)
end

struct ChemPhysics <: AbstractPhysics
    k::Float64
    rho::Float64
    cp::Float64
    chemical_reactions::Vector{ChemicalReaction}
    chemical_vol_source_term::Vector{Float64} #chemical addition
    heat_vol_source_term::Float64 #volumetric heating
end

struct HeatPhysics <: AbstractPhysics
    k::Float64 
    rho::Float64
    cp::Float64
    heat_vol_source_term::Float64 #volumetric heating
end

abstract type AbstractBC end

struct HeatBC <: AbstractBC
    initial_temp::Float64
end

struct ChemBC <: AbstractBC
    initial_mass_fractions::Vector{Float64}
end

struct MultiPhysicsBCs
    chem_bcs::Vector{ChemBC}
    temp_bcs::Vector{HeatBC}
end

struct BoundarySystem
    boundary_map::MultiPhysicsBCs
    free_idxs::Vector{Int}
    dirichlet_idxs::Vector{Int}
end

function get_k_effective(k_a, k_b)
    return 2 * k_a * k_b / (k_a + k_b)
end

function numerical_flux(k_avg, T_L, T_R, area, dist)
    grad_T = (T_R - T_L) / dist
    q = -k_avg * grad_T
    return q * area
end

function diffusion_temp_exchange!(
        du_temp_a, du_temp_b,
        k_a, k_b,
        temp_a, temp_b,
        connection_area, connection_distance,
    )
    
    k_effective = get_k_effective(k_a, k_b)

    F = numerical_flux(k_effective, temp_a, temp_b, connection_area, connection_distance)
    
    du_temp_a -= F
    du_temp_b += F
end
    

function net_reaction_rate(chemical_reaction, molar_concentrations, T, kf_A, kf_Ea, kr_A, kr_Ea)
    R_gas_kj = 0.008314
    kf = (kf_A * exp(-kf_Ea / (R_gas_kj * T)))
    kr = (kr_A * exp(-kr_Ea / (R_gas_kj * T)))

    forward_term = 1.0
    for (i, species_id) in enumerate(chemical_reaction.reactants)
        concentration = molar_concentrations[species_id]
        stoich_coeff = chemical_reaction.reactant_stoich_coeffs[i]
        forward_term *= concentration^stoich_coeff
    end
    
    reverse_term = 1.0
    for (i, species_id) in enumerate(chemical_reaction.products)
        concentration = molar_concentrations[species_id]
        stoich_coeff = chemical_reaction.product_stoich_coeffs[i]
        reverse_term *= concentration^stoich_coeff
    end

    net_reaction_rate = ((kf * forward_term) - (kr * reverse_term))

    return net_reaction_rate
end

function K_gibbs_free(T_ref, T_actual, ΔG_rxn_ref, ΔH_rxn_ref)
    K_ref = exp(-ΔG_rxn_ref / (8.314e-3 * T_ref)) #R is in kJ

    ln_K_ratio = (-ΔH_rxn_ref / 8.314e-3) * (1/T_actual - 1/T_ref)
    
    K_T = K_ref * exp(ln_K_ratio)
    
    return K_T 
end

function react_cell!(
    du_mass_fractions, du_temp, molar_concentrations_cache, net_rates_cache, #mutated vars
    rho, vol, #properties
    species_mass_fractions, cell_temp, #u data
    species_molecular_weights, cell_chemical_reactions_vec #other data
    )

    for (species_id, species_mass_fraction) in enumerate(species_mass_fractions)
        molar_concentrations_cache[species_id] = (rho * species_mass_fraction) / species_molecular_weights[species_id]
    end

    for (reaction_id, reaction) in enumerate(cell_chemical_reactions_vec)
        kf_A = reaction.kf_A
        kf_Ea = reaction.kf_Ea

        #find reverse pre exponential_factor
        K_ref = K_gibbs_free(reaction.K_gibbs_free_ref_temp, cell_temp, reaction.delta_gibbs_free_energy, reaction.heat_of_reaction)

        kr_A = (kf_A / K_ref) * exp(-reaction.heat_of_reaction / (8.314e-3 * cell_temp))

        #find reverse Ea
        kr_Ea = kf_Ea - reaction.heat_of_reaction

        net_rates_cache[reaction_id] = net_reaction_rate(reaction, molar_concentrations_cache, cell_temp, kf_A, kf_Ea, kr_A, kr_Ea)
    end

    for (reaction_id, reaction) in enumerate(cell_chemical_reactions_vec) 
        for (species_id, species_molar_concentration) in enumerate(molar_concentrations_cache)
            change_in_species_molar_concentration = 0.0

            for (reaction_idx, reaction) in enumerate(cell_chemical_reactions_vec)
                stoich = reaction.all_stoich_coeffs[species_id]
                change_in_species_molar_concentration += net_rates_cache[reaction_idx] * stoich

                du_temp[1] += net_rates_cache[reaction_idx] * -reaction.heat_of_reaction * vol
                #add heat of reaction to cell temp, - because -1.0 = exothermic
                #rate is in mol/m^3
                #du.temp[cell_id] is viewed as a 1 dimensional array
            end

            du_mass_fractions[species_id] += (change_in_species_molar_concentration * species_molecular_weights[species_id]) / rho
        end
    end
end


function FVM_iter_f!(
        du, u, p, t,
        cell_volumes, 
        cell_centroids, connection_areas, connection_normals, 
        #cell volumes and cell centroids are accessed at the id of the cell
        connection_distances, unconnected_areas,
        #connection areas, normals, and distances are simply accessed by their location in the 
        #list which corresponds to the respective connection in cell_neighbor_map
        species_molecular_weights,
        cell_props_id_map, bc_sys::BoundarySystem, chem_phys::Vector{ChemPhysics}, #heat_phys::Vector{HeatPhysics}, 
        molar_concentrations_cache, net_rates_cache,
        ax, n_reactions, n_species
    )

    #A_Ea_pairs = eachcol(reshape(p, :, n_reactions))
    #unflattened_p would be [[reaction_1_kf_A, reaction_1_kf_Ea], [reaction_2_kf_A, reaction_2_kf_Ea], etc..] 

    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)
    
    du .= 0.0

    #connections loop
    for (i, (idx_a, idx_b)) in enumerate(cell_neighbor_map)
        prop_a = cell_props_id_map[idx_a]
        prop_b = cell_props_id_map[idx_b]

        k_a = chem_phys[prop_a].k #maybe use heat_phys later
        k_b = chem_phys[prop_b].k

        connection_area = connection_areas[i]
        connection_distance = connection_distances[i]

        diffusion_temp_exchange!(
            du.temp[idx_a], du.temp[idx_b],
            k_a, k_b,
            u.temp[idx_a], u.temp[idx_b],
            connection_area, connection_distance,
        )
    end

    molar_concentrations_cache = zeros(eltype(u.mass_fractions), length(u.mass_fractions[:, 1])) #just using mass fractions for cell 1, this may cause some issues later!
    net_rates_cache = zeros(eltype(u.mass_fractions), length(chem_phys[1].chemical_reactions))

    # Source and Capacity Loop
    for cell_id in bc_sys.free_idxs
        #basic props
        vol = cell_volumes[cell_id]
        props = cell_props_id_map[cell_id]

        k = chem_phys[props].k
        rho = chem_phys[props].rho
        cp  = chem_phys[props].cp
        cell_chemical_reactions_vec = chem_phys[props].chemical_reactions



        #chemical reactions loop
        #species_mass_fractions = u.mass_fractions[:, cell_id]
        cell_temp = u.temp[cell_id]
        species_mass_fractions = view(u.mass_fractions, :, cell_id) #we should maybe use views here, probably does matter that much
        
        react_cell!(
            @view(du.mass_fractions[:, cell_id]), @view(du.temp[cell_id:cell_id]), molar_concentrations_cache, net_rates_cache, #mutated vars
            rho, vol,#properties
            species_mass_fractions, cell_temp, #u data
            species_molecular_weights, cell_chemical_reactions_vec #other data
        )

        # heat source and capacity loop
        S = chem_phys[props].heat_vol_source_term * vol 
        # we should probably create separate containers in heat_phys for both source terms on a per area and per cell basis

        du.temp[cell_id] += S
        
        cap = rho * cp * vol
        du.temp[cell_id] /= cap
    end

    for cell_id in bc_sys.dirichlet_idxs
        du.temp[cell_id] = 0.0
    end
end

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

n_species = length(initial_mass_fractions)
alloc_mass_fraction_vec = zeros(n_species, n_cells) 

u_proto = ComponentArray(mass_fractions = alloc_mass_fraction_vec, temp = zeros(n_cells))

for cell_id in 1:n_cells
    u_proto.mass_fractions[:, cell_id] = chem_bcs[cell_id].initial_mass_fractions
    u_proto.temp[cell_id] = heat_bcs[cell_id].initial_temp
end

cell_props_id_map = Int[]

for cell in CellIterator(grid)
    cell_id = cellid(cell)
    for i in eachindex(cell_sets)
        if cell_id in cell_sets[i].cell_set_idxs
            push!(cell_props_id_map, cell_sets[i].props_id)
        end
    end
end

u_axes = getaxes(u_proto)[1] 
#we use u_proto and u_axes in the function because KrylovJL_GMRES complains when a component array from ComponentArrays.jl is passed in

u0 = Vector(u_proto)

initial_node_coordinates = get_node_coordinates(grid)

cell_neighbor_map, neighbor_map_respective_node_ids = get_neighbor_map(grid)

unconnected_cell_face_map, unconnected_map_respective_node_ids = get_unconnected_map(grid)

nodes_of_cells = get_nodes_of_cells(grid)

reaction_1_kf_A_guess = 10000.0
reaction_1_kf_Ea_guess = 50.0

proto_p_guess = [reaction_1_kf_A_guess, reaction_1_kf_Ea_guess] #27 element Vector{Vector{}

p_guess = reduce(vcat, proto_p_guess) #81 element Vector{} # we do this because Optimization.solve(...) complains about vectors of vectors

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

#= 
NOTE: We might have to declare these as constant to prevent dynamic dispatch in the future but we're leaving this out for now
const cell_neighbor_map = cell_neighbor_map 
const const_neighbor_map_respective_node_ids = neighbor_map_respective_node_ids
=#

n_reactions = length(reaction_physics.chemical_reactions)
n_species = length(initial_mass_fractions)

du0 = u0 .* 0.0

#=
f_closure = (du, u, p) -> FVM_iter_f!(
    du, u, p, 
    cell_volumes, cell_centroids, connection_areas, connection_normals, connection_distances, 
    unconnected_areas, 
    [0.06005, 0.04607, 0.08811, 0.018015],
    cell_props_id_map, bc_sys, chem_phys_vec, u_axes, n_reactions
)

f_closure(u0 .* 0.0, u0, p_guess)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, p_guess), du0, u0, detector)

nl_func = NonlinearFunction(f_closure, jac_prototype = float.(jac_sparsity))

prob = NonlinearProblem(nl_func, u0, p_guess)

println("sol time")

function algebraicmultigrid(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(convert(AbstractMatrix, W)))
    else
        Pl = Plprev
    end
    Pl, nothing
end

function iluzero(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu0(convert(AbstractMatrix, W))
    else
        Pl = Plprev
    end
    Pl, nothing
end

println("timed forward problem")
#@time sol = solve(prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), saveat=t)

@time sol = solve(prob, NonlinearSolve.NewtonRaphson(concrete_jac = true), p=p_guess)
@VSCodeServer.profview sol = solve(prob, NonlinearSolve.NewtonRaphson(concrete_jac = true), p=p_guess)
=#

# Implicit Solving Stuff Below (just add t after p to the FVM_iter_f! function)

f_closure_implicit = (du, u, p, t) -> FVM_iter_f!(
    du, u, p, t,
    cell_volumes, cell_centroids, connection_areas, connection_normals, connection_distances, unconnected_areas, 
    [0.06005, 0.04607, 0.08811, 0.018015],
    cell_props_id_map, bc_sys, chem_phys_vec, 
    molar_concentrations_cache, net_rates_cache,
    u_axes, n_reactions, n_species
)
#just re-add t to the FVM_iter_f! function above to make it compatible with implicit solving

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0, u0, detector)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0, tMax = 0.0, 10000000000.0 
desired_steps = 10
dt = tMax / desired_steps
tspan = (t0, tMax)
t = t0:dt:tMax;

implicit_prob = ODEProblem(ode_func, u0, tspan, p_guess)

function iluzero(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu0(convert(AbstractMatrix, W))
    else
        Pl = Plprev
    end
    Pl, nothing
end

@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
@VSCodeServer.profview sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))

#algebraicmultigrid is better for very large systems (not sure what the cutoff is, though)
#=
function test_predict(θ)
    sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), p=θ, sensealg=InterpolatingAdjoint(autodiff=AutoEnzyme()))
    return Array(sol), sol.t
end
=#

#SteadyStateProblem with DynamicSS: Good if your system might have stability issues or if you want to leverage your existing ODE infrastructure. 
#Also good if you might later want transient optimization.

#SteadyStateProblem with SSRootfind or NonlinearProblem: Better for pure steady-state. Faster when it works. 
#Requires a good initial guess and a well-conditioned Jacobian.

#=
What I have learned after a lot of debugging 
    - We should probably just give up on fixing zygote and ReverseDiff
    - Instead, we should focus on getting Enzyme to work by preventing any dynamic dispatch because it supports array mutation
    - While we try to get that to work, we should just use AutoForwardDiff as the deform box 

=#
