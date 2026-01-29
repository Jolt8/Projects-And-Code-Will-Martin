using Ferrite
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using SciMLSensitivity
using Optimization, OptimizationPolyalgorithms, Zygote
#using Enzyme
using RecursiveArrayTools
using OptimizationOptimJL
using ILUZero
import AlgebraicMultigrid
#import SparseConnectivityTracer, ADTypes
using NonlinearSolve
import Logging
using ComponentArrays
using StaticArrays
using ProfileView

function get_node_coordinates(grid)
    n_nodes = length(grid.nodes)

    node_coordinates = Vector{SVector{3,Float64}}(undef, n_nodes)

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

    cell_faces_node_ids = NTuple{6,NTuple{4,Int64}}[]
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

    cell_neighbor_pairs = Tuple{Int,Int}[]
    neighbor_map_respective_node_ids = NTuple{4,Int}[]

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

    unconnected_cell_face_map = NTuple{2,Int}[]
    unconnected_map_respective_node_ids = NTuple{4,Int}[]

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
    x = a[2] * b[3] - a[3] * b[2]
    y = a[3] * b[1] - a[1] * b[3]
    z = a[1] * b[2] - a[2] * b[1]
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

    unconnected_areas = Vector{SVector{6,T}}(undef, n_unconnected_faces)
    unconnected_normals = Vector{SVector{6,CoordType}}(undef, n_unconnected_faces)

    for (i, (cell_id, face_idx)) in enumerate(unconnected_cell_face_map)
        areas = MVector{6,T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        normals = MVector{6,CoordType}(
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
    reactants::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    reactant_stoich_coeffs::Vector{Int} #(stoich_coeff_for_reactant_1, stoich_coeff_for_reactant_2, etc...)
    products::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    product_stoich_coeffs::Vector{Int} #[stoich_coeff_for_product_1, stoich_coeff_for_product_2, etc...]
    all_stoich_coeffs::Vector{Int} #(stoich_for_chemical_id_1, stoich_for_chemical_id_2, etc...)
    #=
    commenting these out because there are our optimized_parameters

    kf_A # Pre-exponential factor for forward reaction
    kf_Ea # Activation energy for forward reaction
    kr_A # Pre-exponential factor for reverse reaction
    kr_Ea  # Activation energy for reverse reaction
    =#
end

struct ChemPhysics <: AbstractPhysics
    k::Float64
    rho::Float64
    cp::Float64
    chemical_reactions::Vector{ChemicalReaction}
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

function numerical_flux(k_avg, T_L, T_R, area, dist)
    grad_T = (T_R - T_L) / dist
    q = -k_avg * grad_T
    return q * area
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

    ln_K_ratio = (-ΔH_rxn_ref / 8.314e-3) * (1 / T_actual - 1 / T_ref)

    K_T = K_ref * exp(ln_K_ratio)

    return K_T
end

function FVM_iter_f!(
    du, u, p, t, timestamps, temperatures_vec,
    cell_volumes,
    cell_centroids, connection_areas, connection_normals,
    #cell volumes and cell centroids are accessed at the id of the cell
    connection_distances, unconnected_areas,
    #connection areas, normals, and distances are  accessed by their location in the 
    #list which corresponds to the respective connection in cell_neighbor_map
    species_molecular_weights,
    cell_props_id_map, bc_sys::BoundarySystem, chem_phys::Vector{ChemPhysics}, ax, n_reactions, n_species
)

    A_Ea_pairs = eachcol(reshape(p, :, n_reactions))
    #unflattened_p would be [[reaction_1_kf_A, reaction_1_kf_Ea], [reaction_2_kf_A, reaction_2_kf_Ea], etc..] 

    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)

    du .= 0.0

    #=
    for (i, (idx_a, idx_b)) in enumerate(cell_neighbor_map)
        #diffusion here later
    end
    =#

    # Source and Capacity Loop

    time_idx = argmin(abs.(timestamps .- t))

    for cell_id in bc_sys.free_idxs
        vol = cell_volumes[cell_id]
        props = cell_props_id_map[cell_id]

        k = chem_phys[props].k
        rho = chem_phys[props].rho
        cp = chem_phys[props].cp
        cell_chemical_reactions_vec = chem_phys[props].chemical_reactions

        #S = chem_phys[props].source_term * vol 
        # we should probably create separate containers in chem_phys for both source terms on a per area and per cell basis

        species_mass_fractions = view(u.mass_fractions, :, cell_id)

        species_molar_concentrations = [
            (rho * species_mass_fraction) / species_molecular_weights[species_id]
            for (species_id, species_mass_fraction) in enumerate(species_mass_fractions)
        ]
        #this above causes a lot of GC

        for (reaction_id, reaction) in enumerate(cell_chemical_reactions_vec)
            kf_A = A_Ea_pairs[reaction_id][1]
            kf_Ea = A_Ea_pairs[reaction_id][2]

            #find reverse pre exponential_factor
            K_ref = K_gibbs_free(reaction.K_gibbs_free_ref_temp, temperatures_vec[time_idx, cell_id], reaction.delta_gibbs_free_energy, reaction.heat_of_reaction)

            kr_A = (kf_A / K_ref) * exp(-reaction.heat_of_reaction / (8.314e-3 * temperatures_vec[time_idx, cell_id]))

            #find reverse Ea
            kr_Ea = kf_Ea - reaction.heat_of_reaction

            net_rates = [net_reaction_rate(reaction, species_molar_concentrations, temperatures_vec[time_idx, cell_id], kf_A, kf_Ea, kr_A, kr_Ea) for reaction in cell_chemical_reactions_vec]

            for (species_id, species_molar_concentration) in enumerate(species_molar_concentrations)
                change_in_species_molar_concentration = 0.0

                for (reaction_idx, reaction) in enumerate(cell_chemical_reactions_vec)
                    stoich = reaction.all_stoich_coeffs[species_id]
                    change_in_species_molar_concentration += net_rates[reaction_idx] * stoich
                end

                du.mass_fractions[species_id, cell_id] = (change_in_species_molar_concentration * species_molecular_weights[species_id]) / rho
            end
        end
    end
end

grid_dimensions = (1, 1, 1)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

addcellset!(grid, "internal_cells", x -> x != "chicken") # all

# acetic acid, ethanol, water, ethyl acetate
initial_mass_fractions = [0.5, 0.5, 0.0, 0.0]

acetic_acid_ethanol_esterification_rxn = ChemicalReaction(
    -4.0, #Delta H
    -4.0, #Delta Gibbs free at ref temp
    298, #ref temp
    [1, 2], #reactant_ids
    [1, 1], #reactant_stoich_coeffs
    [3, 4], #product_ids
    [1, 1], #product_stoich_coeffs
    [-1, -1, 1, 1] #stoich coefficients -1 = reactant, 1 = product
)

reaction_physics = ChemPhysics(0.6e-3, 1000, 4.184, [acetic_acid_ethanol_esterification_rxn])

chem_phys_vec = [reaction_physics]

internal_cell_set_idxs = Set(getcellset(grid, "internal_cells"))

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

u_proto = ComponentArray(mass_fractions=alloc_mass_fraction_vec)

for cell_id in 1:n_cells
    u_proto.mass_fractions[:, cell_id] = chem_bcs[cell_id].initial_mass_fractions
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

#= 
NOTE: We might have to declare these as constant to prevent dynamic dispatch in the future but we're leaving this out for now
const cell_neighbor_map = cell_neighbor_map 
const const_neighbor_map_respective_node_ids = neighbor_map_respective_node_ids
=#

n_reactions = length(reaction_physics.chemical_reactions)
n_species = length(initial_mass_fractions)

du0 = u0 .* 0.0

#Experimental data retrieval and processing

using XLSX
using Unitful

xf = XLSX.readxlsx("C://Users//wille//Desktop//Projects-And-Code-Will-Martin//Excel Projects//esterification_data_processing_for_julia.xlsx")

typeof(float.(vec(xf["40C"]["A2:A10"])))

struct TrialPreprocessedData
    temperature_timestamps::Vector{Float64}
    temperatures::Vector{Float64}
    moles_timestamps::Vector{Float64}
    moles_respective_temp_idxs::Vector{Int64}
    acetic_acid_moles::Vector{Float64}
    ethanol_moles::Vector{Float64}
    ethyl_acetate_moles::Vector{Float64}
    water_moles::Vector{Float64}
end

T1_temp_end_idx = xf["30C Temps"]["I3"]
T1_temp_timestamps = float.(vec(xf["30C Temps"]["K2:K$T1_temp_end_idx"]))
T1_moles_timestamps = float.(vec(xf["30C"]["A3:A11"]))
T1_moles_respective_temp_idxs = findall(x -> x in T1_moles_timestamps, T1_temp_timestamps)
trial_1_data = TrialPreprocessedData(
    #temperature_timestamps
    float.(vec(xf["30C Temps"]["K2:K$T1_temp_end_idx"])),
    #temperatures
    float.(vec(xf["30C Temps"]["L2:L$T1_temp_end_idx"])),

    #moles_timestamps 
    float.(vec(xf["30C"]["A3:A11"])),
    #moles_respective_temp_idxs
    T1_moles_respective_temp_idxs,
    #acetic_acid_moles
    float.(vec(xf["30C"]["G3:G11"])),
    #ethanol_moles
    float.(vec(xf["30C"]["H3:H11"])),
    #ethyl_acetate_moles
    float.(vec(xf["30C"]["J3:J11"])),
    #water_moles
    float.(vec(xf["30C"]["I3:I11"])),
)


T2_temp_end_idx = xf["40C Temps"]["I3"]
T2_temp_timestamps = float.(vec(xf["40C Temps"]["K2:K$T2_temp_end_idx"]))
T2_moles_timestamps = float.(vec(xf["40C"]["A3:A11"]))
T2_moles_respective_temp_idxs = findall(x -> x in T2_moles_timestamps, T2_temp_timestamps)
trial_2_data = TrialPreprocessedData(
    #temperature_timestamps
    float.(vec(xf["40C Temps"]["K2:K$T2_temp_end_idx"])),
    #temperatures
    float.(vec(xf["40C Temps"]["L2:L$T2_temp_end_idx"])),

    #moles_timestamps 
    float.(vec(xf["40C"]["A3:A11"])),
    #moles_respective_temp_idxs
    T2_moles_respective_temp_idxs,
    #acetic_acid_moles
    float.(vec(xf["40C"]["E3:E11"])),
    #ethanol_moles
    float.(vec(xf["40C"]["F3:F11"])),
    #ethyl_acetate_moles
    float.(vec(xf["40C"]["H3:H11"])),
    #water_moles
    float.(vec(xf["40C"]["G3:G11"])),
)

T3_temp_end_idx = xf["50C Temps"]["I3"]
T3_temp_timestamps = float.(vec(xf["50C Temps"]["K2:K$T3_temp_end_idx"]))
T3_moles_timestamps = float.(vec(xf["50C"]["A3:A14"]))
T3_moles_respective_temp_idxs = findall(x -> x in T3_moles_timestamps, T3_temp_timestamps)
trial_3_data = TrialPreprocessedData(
    #temperature_timestamps
    float.(vec(xf["50C Temps"]["K2:K$T3_temp_end_idx"])),
    #temperatures
    float.(vec(xf["50C Temps"]["L2:L$T3_temp_end_idx"])),

    #moles_timestamps
    float.(vec(xf["50C"]["A3:A14"])),
    #moles_respective_temp_idxs
    T3_moles_respective_temp_idxs,
    #acetic_acid_moles
    float.(vec(xf["50C"]["G3:G14"])),
    #ethanol_moles
    float.(vec(xf["50C"]["H3:H14"])),
    #ethyl_acetate_moles
    float.(vec(xf["50C"]["J3:J14"])),
    #water_moles
    float.(vec(xf["50C"]["I3:I14"])),
)

all_trials_data = [trial_1_data, trial_2_data, trial_3_data]

struct Trial
    temperature_timestamps::Vector{Float64}
    temperatures_vec::Matrix{Float64}

    moles_timestamps::Vector{Float64}
    moles_respective_temp_idxs::Vector{Int64}
    mass_fractions_matrix::Array{Float64,3} #[timestamp, reactant_idx]
end

function get_mass_fraction(species_moles, species_molecular_weight, rho, mixture_volume)
    (((species_moles / mixture_volume) * species_molecular_weight) / rho)
end

trials = Trial[]

for trial in all_trials_data
    temp_timestamps = trial.temperature_timestamps .* u"s"

    temperatures_deg_C = trial.temperatures .* u"°C"

    trial_temperatures = temperatures_deg_C .|> u"K"

    moles_timestamps = trial.moles_timestamps
    moles_respective_temp_idxs = trial.moles_respective_temp_idxs

    acetic_acid_moles = trial.acetic_acid_moles .* u"mol"
    ethanol_moles = trial.ethanol_moles .* u"mol"
    ethyl_acetate_moles = trial.ethyl_acetate_moles .* u"mol"
    water_moles = trial.water_moles .* u"mol"

    mixture_volume = 50u"ml"
    mixture_rho = 1000.0u"kg/m^3"

    acetic_acid_mw = 60.05u"g/mol"
    ethanol_mw = 46.07u"g/mol"
    ethyl_acetate_mw = 88.11u"g/mol"
    water_mw = 18.02u"g/mol"

    acetic_acid_mass_fractions = get_mass_fraction.(acetic_acid_moles, acetic_acid_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
    ethanol_mass_fractions = get_mass_fraction.(ethanol_moles, ethanol_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
    ethyl_acetate_mass_fractions = get_mass_fraction.(ethyl_acetate_moles, ethyl_acetate_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
    water_mass_fractions = get_mass_fraction.(water_moles, water_mw, mixture_rho, mixture_volume) .|> u"kg/kg"

    mass_fractions_matrix = zeros(length(moles_timestamps), n_species, n_cells)
    temperature_vec = zeros(length(temp_timestamps), n_cells)

    for i in eachindex(moles_timestamps)
        mass_fractions_matrix[i, :, :] .= [acetic_acid_mass_fractions[i], ethanol_mass_fractions[i], ethyl_acetate_mass_fractions[i], water_mass_fractions[i]]
    end
    #mass_fractions_matrix is formatted as [time_idx, species_idx, cell_idx]

    for i in eachindex(temp_timestamps)
        temperature_vec[i, :] .= ustrip(trial_temperatures[i])
    end
    #temperature_vec is formatted as [time_idx, cell_idx]

    trial = Trial(ustrip.(temp_timestamps), ustrip.(temperature_vec), ustrip.(moles_timestamps), moles_respective_temp_idxs, mass_fractions_matrix)

    push!(trials, trial)
end

#trials[1]
#trials[2]
#trials[3]

errors = zeros(length(trials))

#=
function update_temp!(integrator)
    integrator.u[4] = (integrator.u[4] * m_initial + m_ice) / m_total
end

update_temp_callback = PresetTimeCallback(300.0, update_temp!)
=#

u_axes = getaxes(u_proto)[1]
#we use u_proto and u_axes in the function because KrylovJL_GMRES complains when a component array from ComponentArrays.jl is passed in

u0 = Vector(u_proto)

reaction_1_kf_A_guess = 1516.363624655201
reaction_1_kf_Ea_guess = 62.94210526315789

proto_p_guess = [reaction_1_kf_A_guess, reaction_1_kf_Ea_guess] #27 element Vector{Vector{}

p_guess = reduce(vcat, proto_p_guess) #81 element Vector{} # we do this because Optimization.solve(...) complains about vectors of vectors

mc_proto_length = n_cells * n_species

f_closure_implicit_pre = (du, u, p, t, temperature_timestamps, temperature_vec) -> FVM_iter_f!(
    du, u, p, t, temperature_timestamps, temperature_vec,
    cell_volumes, cell_centroids, connection_areas, connection_normals, connection_distances,
    unconnected_areas,
    [0.06005, 0.04607, 0.08811, 0.018015],
    cell_props_id_map, bc_sys, chem_phys_vec, u_axes, n_reactions, n_species
)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

function iluzero(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu0(convert(AbstractMatrix, W))
    else
        Pl = Plprev
    end
    Pl, nothing
end

trials_mass_fractions_results = []
trials_timestamps = []

u0_for_jac = trials[1].mass_fractions_matrix[1, :, 1]
du0_for_jac = u0_for_jac .* 0.0

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0_for_jac, u0_for_jac, detector)

for trial in trials
    trial_u0 = trial.mass_fractions_matrix[1, :, 1]
    trial_du0 = trial_u0 .* 0.0

    f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trial.temperature_timestamps, trial.temperatures_vec)

    ode_func = ODEFunction(f_closure_implicit, jac_prototype=float.(jac_sparsity))

    trial_t0, trial_tMax = trial.temperature_timestamps[1], trial.temperature_timestamps[end]
    tspan = (trial_t0, trial_tMax)
    implicit_prob = ODEProblem(ode_func, trial_u0, tspan, p_guess)

    @time sol = solve(implicit_prob, FBDF(linsolve=KrylovJL_GMRES(), precs=iluzero, concrete_jac=true), saveat=trial.moles_timestamps)

    push!(trials_mass_fractions_results, sol.u)
    push!(trials_timestamps, sol.t)
end

println(trials_mass_fractions_results[1])
println(trials_timestamps[1])
println(trials_mass_fractions_results[2])
println(trials_timestamps[2])
println(trials_mass_fractions_results[3])
println(trials_timestamps[3])

using DataFrames
ethanol_mass_fractions = []
acetic_acid_mass_fractions = []
ethyl_acetate_mass_fractions = []
water_mass_fractions = []

for mass_fractions in trials_mass_fractions_results[1]
    push!(ethanol_mass_fractions, mass_fractions[1])
    push!(acetic_acid_mass_fractions, mass_fractions[2])
    push!(ethyl_acetate_mass_fractions, mass_fractions[3])
    push!(water_mass_fractions, mass_fractions[4])
end

df1 = DataFrame(AA=trials_timestamps[1], AB=ethanol_mass_fractions, AC=acetic_acid_mass_fractions, AD=ethyl_acetate_mass_fractions, AE=water_mass_fractions)

ethanol_mass_fractions = []
acetic_acid_mass_fractions = []
ethyl_acetate_mass_fractions = []
water_mass_fractions = []
for mass_fractions in trials_mass_fractions_results[2]
    push!(ethanol_mass_fractions, mass_fractions[1])
    push!(acetic_acid_mass_fractions, mass_fractions[2])
    push!(ethyl_acetate_mass_fractions, mass_fractions[3])
    push!(water_mass_fractions, mass_fractions[4])
end

df2 = DataFrame(AA=trials_timestamps[2], AB=ethanol_mass_fractions, AC=acetic_acid_mass_fractions, AD=ethyl_acetate_mass_fractions, AE=water_mass_fractions)

ethanol_mass_fractions = []
acetic_acid_mass_fractions = []
ethyl_acetate_mass_fractions = []
water_mass_fractions = []
for mass_fractions in trials_mass_fractions_results[3]
    push!(ethanol_mass_fractions, mass_fractions[1])
    push!(acetic_acid_mass_fractions, mass_fractions[2])
    push!(ethyl_acetate_mass_fractions, mass_fractions[3])
    push!(water_mass_fractions, mass_fractions[4])
end
df3 = DataFrame(AA=trials_timestamps[3], AB=ethanol_mass_fractions, AC=acetic_acid_mass_fractions, AD=ethyl_acetate_mass_fractions, AE=water_mass_fractions)

#=
XLSX.writetable("esterification_simulation_results_2.xlsx", 
    "T1" => df1, 
    "T2" => df2,
    "T3" => df3
)
=#

println("done")


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

function loss(θ)
    total_loss = 0.0
    for (i, trial) in enumerate(trials)
        trial_u0_proto = ComponentArray(mass_fractions=trial.mass_fractions_matrix[1, :, :], temp=trial.temperatures_vec[1, :]) #get at [t = 0, all_species, all_cells]
        trial_u0 = Vector(trial_u0_proto)
        tspan = (trial.temperature_timestamps[1], trial.temperature_timestamps[end])
        p_adjusted = [exp(θ[1]), θ[2]]

        f_closure_implicit = (du, u, p, t) -> f_closure_implicit_pre(du, u, p, t, trial.temperature_timestamps, trial.temperatures_vec)

        ode_func = ODEFunction(f_closure_implicit)#, jac_prototype=float.(jac_sparsity))

        prob_trial = ODEProblem(ode_func, trial_u0, tspan, p_adjusted)

        sol = solve(
            prob_trial, Rosenbrock23(), p=p_adjusted,
            sensealg=InterpolatingAdjoint(autodiff=AutoForwardDiff()),
            saveat=trial.temperature_timestamps,
            #callback = current_cb
        )
        #Tsit5 is faster for now, we can do Rosenbrock23()

        if sol.retcode != ReturnCode.Success
            return 1e9
        end

        trial_error = 0.0
        for (moles_timestamp_idx, respective_temp_idxs) in enumerate(trial.moles_respective_temp_idxs)
            pred = sol.u[respective_temp_idxs][1:mc_proto_length]
            obs = vec(trial.mass_fractions_matrix[moles_timestamp_idx, :, :])
            trial_error += sum(abs2, pred .- obs)
        end

        total_loss += trial_error
    end
    return total_loss
end

function loss_for_fixed_Ea(ln_A_val, fixed_Ea)
    return loss([ln_A_val[1], fixed_Ea])
end

#=
trial_u0_proto = ComponentArray(mass_fractions = trials[1].mass_fractions_matrix[1, :, :], temp = trials[1].temperatures_vec[1, :]) #get at t = 0
trial_u0 = Vector(trial_u0_proto)
#prob_trial = remake(implicit_prob, u0=trial_u0, p=[50, 10000], tspan = (trials[1].timestamps[1], trials[1].timestamps[end]))
tspan_for_trial = (trials[1].timestamps[1], trials[1].timestamps[end])
prob_trial = remake(implicit_prob, u0=trial_u0, p=[100, 50], tspan = (trials[1].timestamps[1], trials[1].timestamps[end]))
#prob_trial = ODEProblem(ode_func, trial_u0, tspan_for_trial, [10000, 50])
trials[1].timestamps
#@time sol_1 = solve(prob_trial, Rosenbrock23(), saveat=trials[1].timestamps)
@time sol_1 = solve(prob_trial, Tsit5(), saveat=trials[1].timestamps)
#@time sol_1 = solve(prob_trial, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
=#

@time loss([log(1100), 62]) #remember, kf_A, kf_Ea
loss([log(930.0), 61.599])

using ForwardDiff
using Sparspak
#@time ForwardDiff.gradient(loss, [log(10), 60])

LOSS = [] # Loss accumulator
PRED = [] # prediction accumulator
PARS = [] # parameters accumulator
#=
cb = function (state, l)
    display(l)
    display(state.u)
    append!(LOSS, l)
    append!(PARS, [state.u])
    false
    verbose = true
end
=#

SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true #turn to true to debug EnzymeJVP

Logging.disable_logging(Logging.Warn)  # Disable all warnings
#Logging.disable_logging(Logging.Warn - 1)  # enable all warnings

adtype = Optimization.AutoForwardDiff() #TODO: Get stuff working with Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
#optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype, cons_jac_prototype=jac_sparsity)
#adding jac_sparsity kills performance
#update 2: apparently passing the jac_prototype here is actually uncessary and the only place in which you need to pass it is in the original implicit ODE solve

guess_params = Float64[log(1516.363624655201), 62.94210526315789] #remember, kf_A, kf_Ea
#(Ea = 62.94210526315789, A = 1516.363624655201, Loss = 0.03271151459598594)
loss(guess_params)
#ForwardDiff.gradient(loss, guess_params)
lower_bounds = [log(300.0), 10.0]
upper_bounds = [log(1e12), 1000.0]
#@VSCodeServer.profview for a FlameGraph

#@time ForwardDiff.gradient(loss, [log(1000.0), 50.0])
#@VSCodeServer.profview ForwardDiff.gradient(loss, guess_params)

optprob = Optimization.OptimizationProblem(optf, guess_params, lb=lower_bounds, ub=upper_bounds)

@time res = Optimization.solve(
    optprob,
    #callback=cb,
    OptimizationOptimJL.IpNewton(),
    #LBFGS, BFGS, and Fminbox don't work well here and don't return values that are anywhere close to values defined in literature 
    #IPNewton works kinda fine
    f_abstol=1e-4,
    g_abstol=1e-4,
)



optimized_kf_A = exp(res.u[1])
optimized_kf_Ea = res.u[2]

best_params = res.u

H = ForwardDiff.hessian(θ -> loss(θ), best_params)

n_data_points = sum(length(trial.temperature_timestamps) * n_species for trial in trials)
n_params = length(best_params)
σ2 = loss(best_params) / (n_data_points - n_params)
covariance_matrix = 2 * σ2 * inv(H)

std_errors = sqrt.(diag(covariance_matrix))

println("Optimized A: $(exp(best_params[1])) ± $(std_errors[1])")
println("Optimized Ea: $(best_params[2]) ± $(std_errors[2])")

#=
FINAL OPTIMIZED PARAMETER FOR IA (DO NOT TOUCH)
    - kf_A = 1420.7137027799242 ± 0.04496177269267059
    - Ea = 62.96988585067519 ± 0.11841674095409051
=#

#= 
#this took 200 seconds with an initial guess that was already extremely close to the true parameters, 
#if it takes longer when the values are farther away, this will never converge in time
using OptimizationEvolutionary
@time res = Optimization.solve(
    optprob,
    Evolutionary.CMAES(μ = 40, λ = 100),
)
=#
#=
using OptimizationBBO #for BlackBoxOptim
@time res = Optimization.solve(
    optprob, 
    BBO_adaptive_de_rand_1_bin_radiuslimited(), 
)
=#


kf_A_guess = 1e1

guess_lb = 62.5
guess_ub = 63.1

Ea_range = collect(range(guess_lb, guess_ub, 20))

function loss_for_fixed_Ea(ln_A_val, fixed_Ea)
    return loss([ln_A_val[1], fixed_Ea])
end

function profile_Ea_scan(Ea_range)
    results = []

    for Ea_val in Ea_range
        prob_1D = Optimization.OptimizationProblem(
            Optimization.OptimizationFunction((x, p) -> loss_for_fixed_Ea(x, Ea_val), Optimization.AutoForwardDiff()),
            [log(1000.0)], # Initial guess for A
            lb=[log(100.0)], ub=[log(1e12)]
        )

        sol = Optimization.solve(prob_1D, OptimizationOptimJL.BFGS()) # Brent is optimal for 1D

        push!(results, (Ea=Ea_val, A=exp(sol.u[1]), Loss=sol.objective))
        println("Ea: $Ea_val, Opt_A: $(exp(sol.u[1])), Loss: $(sol.objective)")
    end
    return results
end

results = profile_Ea_scan(Ea_range)

A_range = [results[i].Ea for i in eachindex(results)]
loss_range = [results[i].Loss for i in eachindex(results)]

min_idx = argmin(loss_range)

results[min_idx]
results[min_idx].Ea
results[min_idx].A
results[min_idx].Loss

#FINAL PARAMETERS (DO NOT TOUCH)
#61.14478947368421
#650.6316641398093
#0.03102647047431776

#using Plots
plot(Ea_range, A_range)
plot(Ea_range, loss_range)

new_guess_params = [exp(res.u[1]), res.u[2]]
#kf_A is in units of m^3 / (mol * s)

#= Note to future me on doing this
    - Basically, the reason why the Ea and A never really match up with experimental data is because there's too much noise in our titrations
    - the noise in our titrations is bad because the difference in temperatures between the two trials (ΔT ~ 11 K) is too small.
    - This is an issue because we need properly distinguished trials at different temperatures to prevent the solver from increasing Ea or A and in response, decreasing A or Ea
    - Basically, this noise makes the solver have a bunch of possible Ea and A pairs that work to minimize the loss function
    - To decouple Ea and A to fix this, increase the ΔT between trials and reduce noise 
=#

best_params = res.u
H = ForwardDiff.hessian(θ -> loss(θ), best_params)

n_data_points = sum(length(trial.temperature_timestamps) * n_species for trial in trials)
n_params = length(best_params)
σ2 = loss(best_params) / (n_data_points - n_params)
covariance_matrix = 2 * σ2 * inv(H)

#=
n_obs = sum(length(t.timestamps) for t in trials) * n_species
mse = loss(best_params) / (n_obs - length(H))
mse_cov_mat = 2 * mse * inv(H)
=#

std_errors = sqrt.(diag(covariance_matrix))

println("Optimized A: $(exp(best_params[1])) ± $(std_errors[1])")
println("Optimized Ea: $(best_params[2]) ± $(std_errors[2])")

#=
What I have learned after a lot of debugging 
    - We should probably just give up on fixing zygote and ReverseDiff
    - Instead, we should focus on getting Enzyme to work by preventing any dynamic dispatch because it supports array mutation
    - While we try to get that to work, we should just use AutoForwardDiff as the deform box 

=#
