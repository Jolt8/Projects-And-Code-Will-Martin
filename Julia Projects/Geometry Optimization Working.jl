using Ferrite
using DifferentialEquations
using LinearAlgebra
using Dates
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

struct CellData{T, V}
    volume::T
    centroid::V
end

struct Connection{T, V}
    area::T
    normal::V
    distance::T
end

abstract type AbstractBC end

struct HeatBC <: AbstractBC
    type::Symbol 
    initial::Float64
end

abstract type AbstractPhysics end

struct HeatPhysics <: AbstractPhysics
    k::Float64 
    rho::Float64
    cp::Float64
    source_term::Float64 #volumetric heating
end

struct FVMMesh{C, K}
    mesh_cells::C
    connections::K
end

struct MultiPhysicsBCs
    temp1::Vector{HeatBC}
    temp2::Vector{HeatBC} 
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

function source(source_term, vol)
    return source_term * vol
end

function capacity(rho, cp, vol)
    return rho * cp * vol
end

function get_k_effective(k_a, k_b)
    return 2 * k_a * k_b / (k_a + k_b)
end

function build_fvm_geometry(grid, cell_sets)
    n_cells = getncells(grid)

    ref_shape = getrefshape(grid.cells[1])
    poly_interp = Lagrange{ref_shape, 1}() #1 = linear elements 2 = quadratic/curved edges
    facet_qr = FacetQuadratureRule{ref_shape}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation
    facet_values = FacetValues(facet_qr, poly_interp)
    cell_qr = QuadratureRule{ref_shape}(2)
    cell_values = CellValues(cell_qr, poly_interp)

    cells_data = Vector{CellData}(undef, n_cells)
    
    for cell in CellIterator(grid) #You have to use a CellIterator Here
        cell_id = cellid(cell)
        Ferrite.reinit!(cell_values, cell) #not using Ferrite. before reinit breaks this! (DO NOT USE reinit! ON ITS OWN)
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values))
        coords = getcoordinates(cell)
        cent = sum(coords) / length(coords)

        cells_data[cell_id] = CellData(vol, cent)
    end

    top = ExclusiveTopology(grid)
    connections = Vector{Connection}()
    for i in 1:n_cells
        for face_idx in 1:nfacets(grid.cells[i])
            neighbor_info = top.face_face_neighbor[i, face_idx]

            if !isempty(neighbor_info)
                neighbor_idx = collect(neighbor_info[1].idx)[1]
                
                if i < neighbor_idx
                    coords = getcoordinates(grid, i)
                    Ferrite.reinit!(facet_values, coords, face_idx) #not using Ferrite. before reinit breaks this! (DO NOT USE reinit! ON ITS OWN)

                    area = sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values))

                    n_ref = getnormal(facet_values, 1) 

                    centroid_a = cells_data[i].centroid
                    centroid_b = cells_data[neighbor_idx].centroid
                    dist = norm(centroid_b - centroid_a)
                    #TODO: fix for non-orthogonal meshes

                    push!(connections, Connection(area, n_ref, dist))
                end
            end
        end
    end
    
    return FVMMesh(cells_data, connections)
end

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

function get_neighbor_map(grid)
    nodes_of_cells = get_nodes_of_cells(grid)
    #returns:
    #cell - neighbor pairs (ex. (1, 2) (cell_1 is connected with cell_2))
    #respective node idxs of cell and neighbor (ex. (2, 8, 44, 38))

    cell_neighbor_pairs = Tuple{Int, Int}[]
    respective_node_ids = NTuple{4, Int}[]

    #= OLD, BUT MORE READABLE
    cell_neighbor_pairs = Vector[]
    respective_node_ids = Vector[]
    =#

    top = ExclusiveTopology(grid)
    connections = Vector{Connection}()
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
                    push!(respective_node_ids, ntuple(i -> respective_nodes[i], 4))

                    #= OLD BUT MORE READABLE
                    push!(cell_neighbor_pairs, [cell_id, neighbor_id])
                    push!(respective_node_ids, respective_nodes)
                    =#
                end
            end
        end
    end
    return cell_neighbor_pairs, respective_node_ids
end

function precompute_ffd_data(db_coordinates, db_nodes_of_cells, node_coordinates)
    n_nodes = size(node_coordinates, 1)
    
    cell_map = zeros(Int, n_nodes)
    intrinsic_coordinates = Vector{SVector{3, Float64}}(undef, n_nodes)
    
    for node_id in 1:n_nodes
        node_x = node_coordinates[node_id][1]
        node_y = node_coordinates[node_id][2]
        node_z = node_coordinates[node_id][3]
        
        for db_cell_id in eachindex(db_nodes_of_cells)
            curr_cell_nodes = db_nodes_of_cells[db_cell_id]
            
            cell_x_coords = [db_coordinates[n][1] for n in curr_cell_nodes]
            cell_y_coords = [db_coordinates[n][2] for n in curr_cell_nodes]
            cell_z_coords = [db_coordinates[n][3] for n in curr_cell_nodes]
            
            x_min, x_max = extrema(cell_x_coords)
            y_min, y_max = extrema(cell_y_coords)
            z_min, z_max = extrema(cell_z_coords)
            
            eps = 1e-10
            in_x = (x_min - eps) <= node_x <= (x_max + eps)
            in_y = (y_min - eps) <= node_y <= (y_max + eps)
            in_z = (z_min - eps) <= node_z <= (z_max + eps)
            
            if in_x && in_y && in_z
                cell_map[node_id] = db_cell_id
                
                center_x = (x_min + x_max) / 2.0
                center_y = (y_min + y_max) / 2.0
                center_z = (z_min + z_max) / 2.0
                
                half_width_x = (x_max - x_min) / 2.0
                half_width_y = (y_max - y_min) / 2.0
                half_width_z = (z_max - z_min) / 2.0
                
                #avoid division by zero
                if half_width_x > 1e-12
                    xi = (node_x - center_x) / half_width_x
                else
                    xi = 0.0
                end

                if half_width_y > 1e-12
                    eta = (node_y - center_y) / half_width_y
                else
                    eta = 0.0
                end

                if half_width_z > 1e-12
                    zeta = (node_z - center_z) / half_width_z
                else 
                    zeta = 0.0
                end

                intrinsic_coordinates[node_id] = [xi, eta, zeta]
                break
            end
        end
        
        # Warning if node wasn't found in any cell
        if cell_map[node_id] == 0
            @warn "Node $node_id not found in any deform box cell"
        end
    end
    
    return cell_map, intrinsic_coordinates
end

function calculate_single_node_ffd(node_idx, p_control_points, db_nodes_of_cells, cell_map, local_coords)
    # Get the deform box cell this node belongs to
    cell_id = cell_map[node_idx]
    
    # Get local coords
    local_coord = local_coords[node_idx]
    xi = local_coord[1]
    eta = local_coord[2]
    zeta = local_coord[3]
    
    # Get control point indices
    control_point_indices = db_nodes_of_cells[cell_id]
    
    # Reference signs (Ordering must match Hex8 node ordering)
    ref_signs = (
        SVector(-1.0, -1.0, -1.0), SVector( 1.0, -1.0, -1.0),
        SVector( 1.0,  1.0, -1.0), SVector(-1.0,  1.0, -1.0),
        SVector(-1.0, -1.0,  1.0), SVector( 1.0, -1.0,  1.0),
        SVector( 1.0,  1.0,  1.0), SVector(-1.0,  1.0,  1.0)
    )
    
    val_x, val_y, val_z = 0.0, 0.0, 0.0
    
    for corner in 1:8
        signs = ref_signs[corner]

        weight = 0.125 * (1.0 + signs[1] * xi) * 
                         (1.0 + signs[2] * eta) * 
                         (1.0 + signs[3] * zeta)
        
        cp_idx = control_point_indices[corner]
        cp = p_control_points[cp_idx]
        
        val_x += weight * cp[1]
        val_y += weight * cp[2]
        val_z += weight * cp[3]
    end
    
    return SVector(val_x, val_y, val_z)
end

#=
function apply_ffd_motion(p_control_points, db_nodes_of_cells, cell_map, local_coords)
    n_fine_nodes = size(local_coords, 1)

    new_positions_vec = map(1:n_fine_nodes) do i
        calculate_single_node_ffd(i, p_control_points, db_nodes_of_cells, cell_map, local_coords)
    end
    
    return new_positions_vec
end
=#

#no map version

function apply_ffd_motion(p_control_points, db_nodes_of_cells, cell_map, local_coords)
    n_fine_nodes = size(local_coords, 1)

    T = eltype(eltype(p_control_points))

    new_positions_vec = Vector{SVector{3, T}}(undef, n_fine_nodes)

    for i in 1:n_fine_nodes
        new_positions_vec[i] = calculate_single_node_ffd(i, p_control_points, db_nodes_of_cells, cell_map, local_coords)
    end
    
    return new_positions_vec
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

#=
function rebuild_fvm_geometry(cell_neighbor_map, respective_node_ids, node_coordinates, nodes_of_cells)
    #The cells_data and connections map are still causing GC
    cells_data = map(eachindex(nodes_of_cells)) do cell_id
        cell_nodes = nodes_of_cells[cell_id]
        
        #this is ugly, but it avoids allocation
        p = ntuple(8) do i
            @inbounds node_coordinates[cell_nodes[i]]
        end

        #= more readable:
        p = Vector{SVector{3}}(undef, n_nodes_in_cell)

        for i in eachindex(cell_nodes)
            p[i] = node_coordinates[cell_nodes[i]]
        end
        =#
        
        vol = calculate_hex_volume(p)
        
        cent = sum(p) / 8.0

        CellData(vol, cent)
    end

    connections = map(enumerate(cell_neighbor_map)) do (i, (cell_id, neighbor_id))
        face_node_indices = respective_node_ids[i]
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
        dist = norm(cells_data[cell_id].centroid - cells_data[neighbor_id].centroid)

        #Even if we comment out everything from connections = map()... to this comment, we still get GC for some reason

        Connection(total_area, cell_normal, dist)
    end
    
    return FVMMesh(cells_data, connections)
end
=#


function rebuild_fvm_geometry(cell_neighbor_map, respective_node_ids, node_coordinates, nodes_of_cells)
    #The cells_data and connections map are still causing GC
    CoordType = eltype(node_coordinates)
    T = eltype(CoordType)

    n_cells = length(nodes_of_cells)

    cells_data = Vector{CellData{T, CoordType}}(undef, n_cells)

    for cell_id in eachindex(nodes_of_cells)
        cell_nodes = nodes_of_cells[cell_id]

        n_nodes = length(cell_nodes)

        p = ntuple(8) do i
            @inbounds node_coordinates[cell_nodes[i]]
        end

        vol = calculate_hex_volume(p)
        
        cent = sum(p) / T(n_nodes)
        
        cells_data[cell_id] = CellData(vol, cent)
    end

    n_connections = length(cell_neighbor_map)

    connections = Vector{Connection{T, CoordType}}(undef, n_connections)

    for (i, (cell_id, neighbor_id)) in enumerate(cell_neighbor_map)
        face_node_indices = respective_node_ids[i]
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
        dist = norm(cells_data[cell_id].centroid - cells_data[neighbor_id].centroid)

        #Even if we comment out everything from connections = map()... to this comment, we still get GC for some reason

        connections[i] = Connection(total_area, cell_normal, dist)
    end
    
    return FVMMesh(cells_data, connections)
end

function FVM_iter_f!(
        du, u, p,
        cell_neighbor_map, respective_node_ids, nodes_of_cells, db_nodes_of_cells,
        cell_map, intrinsic_coordinates,
        cell_mat_id_map, bc_sys::BoundarySystem, heat_phys::Vector{HeatPhysics}, ax, db_grid_n_nodes
    )
    #you could add another vector for fluid_phys if needed like fluid_phys::Vector{FluidPhysics}

    unflattened_p = eachcol(reshape(p, 3, db_grid_n_nodes))

    new_node_coordinates = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)

    geo = rebuild_fvm_geometry(cell_neighbor_map, respective_node_ids, new_node_coordinates, nodes_of_cells)
    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)
    
    du .= 0.0

    for (i, (idx_a, idx_b)) in enumerate(cell_neighbor_map)
        conn = geo.connections[i]

        mat_a = cell_mat_id_map[idx_a]
        mat_b = cell_mat_id_map[idx_b]

        k_a = heat_phys[mat_a].k
        k_b = heat_phys[mat_b].k

        k_effective = get_k_effective(k_a, k_b)

        #temp 1 
        u_a = u.temp1[idx_a]
        u_b = u.temp1[idx_b]

        F_1 = numerical_flux(k_effective, u_a, u_b, conn.area, conn.distance)
        
        du.temp1[idx_a] -= F_1
        du.temp1[idx_b] += F_1

        #temp 2 - while this variable seems useless, it's just to test multiple variables
        u_a_2 = u.temp2[idx_a]
        u_b_2 = u.temp2[idx_b]

        F_2 = numerical_flux(k_effective, u_a_2, u_b_2, conn.area, conn.distance)

        du.temp2[idx_a] -= F_2
        du.temp2[idx_b] += F_2
    end

    # Source and Capacity Loop
    for i in bc_sys.free_idxs
        vol = geo.mesh_cells[i].volume
        mat = cell_mat_id_map[i]
        
        rho = heat_phys[mat].rho
        cp  = heat_phys[mat].cp

        S = heat_phys[mat].source_term * vol 
        # we should probably create separate containers in heat_phys for both source terms on a per area and per cell basis

        du.temp1[i] += S
        du.temp2[i] += S
        #=
        h_coeff = 50.0 # Convection coefficient

        T_ambient = 298.0
        area = geo.mesh_cells[i].volume^(2/3) # Approximation of surface area for a cell, or calculate properly

        # Heat loss to environment
        q_convection = h_coeff * area * (u.temp1[i] - T_ambient)

        # Subtract this from the residual
        du.temp1[i] -= q_convection
        =#
        
        cap = capacity(rho, cp, vol)
        du.temp1[i] /= cap
        du.temp2[i] /= cap
    end

    for i in bc_sys.dirichlet_idxs
        du.temp1[i] = 0.0
        du.temp2[i] = 0.0
    end
end

grid_dimensions = (10, 5, 5) 
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

cell_half_dist = (right[1] / grid_dimensions[1]) / 2 
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))

heated_copper_physics = HeatPhysics(401.0, 8960.0, 385.0, 0) #note that if we created a new struct for copper_physics performance would die
steel_physics = HeatPhysics(30.0, 8000.0, 460.0, 0)

heat_phys_vec = [heated_copper_physics, steel_physics]

copper_cell_set_idxs = Set(getcellset(grid, "copper"))
steel_cell_set_idxs = Set(getcellset(grid, "steel"))

struct CellSet
    mat_id::Int
    cell_set_idxs::Set{Int}
end

copper_set = CellSet(1, copper_cell_set_idxs)
steel_set = CellSet(2, steel_cell_set_idxs)

cell_sets = [copper_set, steel_set]

free_idxs = Int[]
dirichlet_idxs = Int[]

function my_bc_mapper(cell_id)
    if cell_id in copper_cell_set_idxs
        bc_type_a = HeatBC(:Neumann, 500.0) #use :Dirichlet to fix temperature to initial in HeatBC
        bc_type_b = HeatBC(:Neumann, 700.0)
        push!(free_idxs, cell_id)
        return [bc_type_a, bc_type_b] 
    elseif cell_id in steel_cell_set_idxs
        bc_type_a = HeatBC(:Neumann, 300.0)
        bc_type_b = HeatBC(:Neumann, 400.0)
        push!(free_idxs, cell_id)
        return [bc_type_a, bc_type_b]
    end
end

temp1_bcs = HeatBC[]
temp2_bcs = HeatBC[]

n_cells = length(grid.cells)
n_vars = 2

for cell_id in 1:n_cells
    bcs = my_bc_mapper(cell_id) #returns vector [BC1, BC2]
    push!(temp1_bcs, bcs[1])
    push!(temp2_bcs, bcs[2])
end

boundary_map = MultiPhysicsBCs(temp1_bcs, temp2_bcs)

bc_sys = BoundarySystem(boundary_map, free_idxs, dirichlet_idxs)

initial_fvm_geo = build_fvm_geometry(grid, cell_sets)

db_grid_dimensions = (5, 5, 5)
db_left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
db_right = Ferrite.Vec{3}((1.1, 1.1, 1.1)) #has to be larger than the original grid
db_grid = generate_grid(Hexahedron, db_grid_dimensions, db_left, db_right)

u_proto = ComponentArray(temp1 = zeros(n_cells), temp2 = zeros(n_cells))

for cell_id in 1:n_cells
    u_proto.temp1[cell_id] = bc_sys.boundary_map.temp1[cell_id].initial
    u_proto.temp2[cell_id] = bc_sys.boundary_map.temp2[cell_id].initial
end

cell_mat_id_map = Int[]

for cell in CellIterator(grid)
    cell_id = cellid(cell)
    for i in eachindex(cell_sets)
        if cell_id in cell_sets[i].cell_set_idxs
            push!(cell_mat_id_map, cell_sets[i].mat_id)
        end
    end
end

u_axes = getaxes(u_proto)[1] 
#we use u_proto and u_axes in the function because KrylovJL_GMRES complains when a component array from ComponentArrays.jl is passed in

u0 = Vector(u_proto)

initial_node_coordinates = get_node_coordinates(grid)

cell_neighbor_map, respective_node_ids = get_neighbor_map(grid)

nodes_of_cells = get_nodes_of_cells(grid)

proto_p_guess = get_node_coordinates(db_grid) #27 element Vector{Vector{}

p_guess = reduce(vcat, proto_p_guess) #81 element Vector{} # we do this because Optimization.solve(...) complains about vectors of vectors

db_nodes_of_cells = get_nodes_of_cells(db_grid)

cell_map, intrinsic_coordinates = precompute_ffd_data(proto_p_guess, db_nodes_of_cells, initial_node_coordinates)

db_grid_n_nodes = length(db_grid.nodes)

unflattened_p = eachcol(reshape(p_guess, 3, db_grid_n_nodes))

new_node_coordinates = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)

const const_cell_neighbor_map = cell_neighbor_map
const const_respective_node_ids = respective_node_ids
const const_nodes_of_cells = nodes_of_cells
const const_db_nodes_of_cells = db_nodes_of_cells

const const_cell_map = cell_map
const const_intrinsic_coordinates = intrinsic_coordinates

const const_cell_mat_id_map = cell_mat_id_map
const const_bc_sys = bc_sys
const const_heat_phys_vec = heat_phys_vec
const const_u_axes = u_axes
const const_db_grid_n_nodes = db_grid_n_nodes

f_closure = (du, u, p) -> FVM_iter_f!(
    du, u, p, 
    const_cell_neighbor_map, const_respective_node_ids, const_nodes_of_cells, const_db_nodes_of_cells,
    const_cell_map, const_intrinsic_coordinates, 
    const_cell_mat_id_map, const_bc_sys, const_heat_phys_vec, const_u_axes, const_db_grid_n_nodes
)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

du0 = u0 .* 0.0

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

no = 1

# Implicit Solving Stuff Below (just add t after p to the FVM_iter_f! function)
#=
f_closure_implicit = (du, u, p, t) -> FVM_iter_f!(
    du, u, p, t,
    cell_neighbor_map, respective_node_ids, cell_sets, nodes_of_cells, db_nodes_of_cells,
    cell_map, intrinsic_coordinates,
    cell_mat_id_map, bc_sys, heat_phys_vec, u_axes, db_grid_n_nodes
)
#just re-add t to the FVM_iter_f! function above to make it compatible with implicit solving

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0, u0, detector)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0, tMax = 0.0, 1000.0 
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
=#
#algebraicmultigrid is better for very large systems (not sure what the cutoff is, though)


function test_predict(θ)
    sol = solve(prob, NonlinearSolve.NewtonRaphson(), p=θ, sensealg=SteadyStateAdjoint(autodiff=AutoEnzyme()), verbose = true, alias=SciMLBase.NonlinearAliasSpecifier(alias_u0=true))
    return Array(sol), SciMLBase.successful_retcode(sol)
end

#SteadyStateProblem with DynamicSS: Good if your system might have stability issues or if you want to leverage your existing ODE infrastructure. 
#Also good if you might later want transient optimization.

#SteadyStateProblem with SSRootfind or NonlinearProblem: Better for pure steady-state. Faster when it works. 
#Requires a good initial guess and a well-conditioned Jacobian.

@time test_pred = test_predict(p_guess)[1]

desired_temp = 500.0

grid_idx_closest_to_desired_temp = argmin(abs.(test_pred .- desired_temp))

function physics_loss(pred_array)
    final_temp = pred_array[grid_idx_closest_to_desired_temp]
    physics_penalty = abs2(final_temp - desired_temp)
    return physics_penalty
end

initial_physics_loss = physics_loss(test_pred)

LOSS = [] # Loss accumulator
PRED = [] # prediction accumulator
PARS = [] # parameters accumulator

cb = function (state, l)
    display(l)
    display(state.u)
    pred = predict(state.u)
    append!(PRED, [pred])
    append!(LOSS, l)
    append!(PARS, [state.u])
    false
    verbose = true
end

SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = true #turn to true to debug EnzymeJVP

lower_bounds = zeros(n_cells)
upper_bounds = ones(n_cells)

reuse_prev_optimized_parameters = false

if reuse_prev_optimized_parameters 
    initial_parameters = new_guess_params
else 
    initial_parameters = p_guess
end

new_guess_params = 0

Logging.disable_logging(Logging.Warn)  # Disable all warnings
#Logging.disable_logging(Logging.Warn - 1)  # enable all warnings

function predict(θ)
    #not sure if adding the linsolve = ... even does anything here, I've tested it and it doesn't seem to make that much of a difference
    #sol = solve(prob, NonlinearSolve.NewtonRaphson(), p=θ, sensealg=SteadyStateAdjoint(diff_type=AutoEnzyme(), linsolve = KrylovJL_GMRES()), verbose = false, alias=SciMLBase.NonlinearAliasSpecifier(alias_u0=true))
    sol = solve(prob, NonlinearSolve.NewtonRaphson(), p=θ, sensealg=SteadyStateAdjoint(diff_type = AutoEnzyme()), verbose = true, alias=SciMLBase.NonlinearAliasSpecifier(alias_u0=true))
    return sol#, SciMLBase.successful_retcode(sol)
end

function get_total_volume(p::Vector{Float64})
    unflattened_p = eachcol(reshape(p, 3, db_grid_n_nodes))
    new_node_coordinates = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)
    
    # We only need the cells_data part of rebuild_fvm_geometry
    total_vol = 0.0
    for cell_id in 1:n_cells
        cell_nodes = nodes_of_cells[cell_id]
        p_vec = ntuple(Val(8)) do i
            new_node_coordinates[cell_nodes[i]]
        end
        total_vol += calculate_hex_volume(p_vec)
    end
    return total_vol
end

function loss(θ)
    #pred_u, is_successful = predict(θ) #if we ever need to return a high value for crashing the solver
    pred_u = predict(θ)
    
    scaled_physics_loss = physics_loss(pred_u) / initial_physics_loss 

    #current_vol = get_total_volume(θ) # first arg is ignored in our helper
    #println(current_vol)
    #vol_loss = abs2(current_vol - initial_volume)

    #println("scaled_physics_loss, ", scaled_physics_loss)

    return (1.0 * scaled_physics_loss) #+ (vol_loss)
end

new_guess_params = initial_parameters

adtype = Optimization.AutoForwardDiff() #TODO: Get stuff working with Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
#optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype, cons_jac_prototype=jac_sparsity)
#adding jac_sparsity kills performance
#update 2: apparently passing the jac_prototype here is actually uncessary and the only place in which you need to pass it is in the original implicit ODE solve

optprob = Optimization.OptimizationProblem(optf, new_guess_params)

#@VSCodeServer.profview for a FlameGraph

@VSCodeServer.profview res = Optimization.solve(
    optprob, OptimizationOptimJL.BFGS(), 
    callback=cb,
    f_abstol = 1e-6, 
    g_abstol = 1e-6,
    maxiters = 100
)

new_guess_params = res.u

res.u == p_guess

db_n_nodes = length(db_grid.nodes)

optimized_parameters = eachcol(reshape(res.u, 3, db_n_nodes))

println(optimized_parameters)

loss(p_guess)

loss(res.u)

test_change = copy(res.u)

test_change[1:54] .+= 1.0

test_change[70:end] .+= 2.0

test_change

loss(test_change)

get_total_volume(test_change)

unflattened_p = eachcol(reshape(test_change, 3, db_grid_n_nodes))

initial_node_coordinates

new_node_coordinates = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)

geo = rebuild_fvm_geometry(cell_neighbor_map, respective_node_ids, new_node_coordinates, nodes_of_cells)

geo.mesh_cells[1].volume

geo.connections[1].distance



#=
TODO list
    - Get rid of geo in favour of just passing back a bunch of lists in rebuild_fvm_geometry
    - Use ComponentArrays to handle instances of metal wall/fluid cell interfaces as they are techically connected, but shouldn't be considered connected
    - in rebuild_fvm_geometry, 
    - Pass the lists found in rebuild_fvm_geometry after implementing the above to avoid having to create new vectors every iteration and triggering GC
    - Implement FEM to calculate wall strain for the microchannel SMR
    - Deforming the mesh based on the FEM solve is not worth implementing at least for right now
=#




#=
What I have learned after a lot of debugging 
    - We should probably just give up o fixing zygote and ReverseDiff
    - Instead, we should focus on getting Enzyme to work by preventing any dynamic dispatch because it supports array mutation
    - While we try to get that to work, we should just use AutoForwardDiff as the deform box 

=#

#=
#Nothing below returns NaN and all of them return gradients that make sense 
Zygote.gradient(p -> begin
    unflattened_p = eachcol(reshape(p, 3, db_grid_n_nodes))
    new_coords = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)
    rebuild_fvm_geometry(cell_neighbor_map, respective_node_ids, new_coords, nodes_of_cells).mesh_cells[1].volume
end, p_guess)[1]

Zygote.gradient(p -> begin
    unflattened_p = eachcol(reshape(p, 3, db_grid_n_nodes))
    new_coords = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)
    rebuild_fvm_geometry(cell_neighbor_map, respective_node_ids, new_coords, nodes_of_cells).mesh_cells[1].centroid[1]
end, p_guess)[1]

Zygote.gradient(p -> begin
    unflattened_p = eachcol(reshape(p, 3, db_grid_n_nodes))
    new_coords = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)
    rebuild_fvm_geometry(cell_neighbor_map, respective_node_ids, new_coords, nodes_of_cells).connections[1].area
end, p_guess)[1] 
#this one returns a 2 gradients that are close to 1.0e-17

Zygote.gradient(p -> begin
    unflattened_p = eachcol(reshape(p, 3, db_grid_n_nodes))
    new_coords = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)
    rebuild_fvm_geometry(cell_neighbor_map, respective_node_ids, new_coords, nodes_of_cells).connections[1].normal[3]
end, p_guess)[1] 
#this one returns a 4 gradients that are close to 1.0e-17

Zygote.gradient(p -> begin
    unflattened_p = eachcol(reshape(p, 3, db_grid_n_nodes))
    new_coords = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)
    rebuild_fvm_geometry(cell_neighbor_map, respective_node_ids, new_coords, nodes_of_cells).connections[1].distance
end, p_guess)[1]

Zygote.gradient(p -> begin
    unflattened_p = eachcol(reshape(p, 3, db_grid_n_nodes))
    new_coords = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)
    sum(sum.(new_coords))
end, p_guess)[1] #this one works just fine

Zygote.gradient(p -> get_total_volume(nothing, p), p_guess) #this one also works just fine

Zygote.gradient(p -> loss(p), p_guess)[1] #this is the only one that returns NaN

Zygote.gradient(p -> sum(predict(p)), p_guess)
=#
#=
function loss_enzyme(θ)
    pred_u = predict(θ)
    scaled_physics_loss = physics_loss(pred_u) / initial_physics_loss 
    #current_vol = get_total_volume(θ)
    #vol_loss = abs2(current_vol - initial_volume)
    return (1.0 * scaled_physics_loss) #+ vol_loss
end

# Use Enzyme directly instead of Zygote
d_p = Enzyme.make_zero(p_guess)
try
    Enzyme.autodiff(Reverse, loss_enzyme, Active, Duplicated(p_guess, d_p))
    println("It Worked!")
catch err
    code_typed(err)
end
=#
#Tried doing this with an explicitly solved problem and I got the same error message, code_typed doesn't even do anything 

