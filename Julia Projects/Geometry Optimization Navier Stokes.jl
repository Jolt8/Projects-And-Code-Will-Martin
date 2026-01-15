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

grid_dimensions = (3, 3, 3) 
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

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

#=
NOTE TO FUTURE SELF ON PURE MUTATION:
    - Check Geometry Optimization Notes.txt or Geometry Optimization Only Mutation.jl for a full write-up 
    - TLDR: don't pass in a geometry cache or node_coordinates cache into FVM_iter_f! to avoid GC in creating new vectors in rebuild_fvm_geometry
=#


abstract type AbstractBC end

struct VelBC <: AbstractBC
    type::Symbol 
    initial::Float64
end

struct PressureBC <: AbstractBC
    type::Symbol 
    initial::Float64
end

abstract type AbstractPhysics end

struct FluidPhysics <: AbstractPhysics
    k::Float64 
    rho::Float64
    mu::Float64
    cp::Float64
end

struct MultiPhysicsBCs
    vel_x::Vector{VelBC}
    vel_y::Vector{VelBC}
    vel_z::Vector{VelBC}
    pressure::Vector{PressureBC}
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

function get_prop_harmonic_mean(prop_a, prop_b)
    return 2 * prop_a * prop_b / (prop_a + prop_b)
end

function upwind(u_left, u_right, mass_flow_rate)
    if mass_flow_rate > 0
        return u_left
    else 
        return u_right
    end
end

[Float64[0.0, 0.0, 0.0] for _ in 1:10]

[Real[0.0, 0.0, 0.0] for _ in 1:10]

function FVM_iter_f!(
        du, u, p,
        cell_neighbor_map, neighbor_map_respective_node_ids, 
        unconnected_cell_face_map, unconnected_map_respective_node_ids,
        nodes_of_cells, db_nodes_of_cells,
        cell_map, intrinsic_coordinates,
        cell_mat_id_map, bc_sys::BoundarySystem, fluid_phys::Vector{FluidPhysics}, ax, db_grid_n_nodes
    )
    #you could add another vector for fluid_phys if needed like fluid_phys::Vector{FluidPhysics}

    unflattened_p = eachcol(reshape(p, 3, db_grid_n_nodes))

    new_node_coordinates = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)

    cell_volumes, 
    cell_centroids, #cell volumes and cell centroids are accessed at the id of the cell
    connection_areas, 
    connection_normals, 
    connection_distances, #connection areas, normals, and distances are simply accessed by their location in the list which corresponds to the respective connection in cell_neighbor_map
    unconnected_areas,
    unconnected_normals = rebuild_fvm_geometry(
        cell_neighbor_map, neighbor_map_respective_node_ids, 
        unconnected_cell_face_map, unconnected_map_respective_node_ids,
        new_node_coordinates, nodes_of_cells
    )

    #there's gotta be a better way to format this shit above
    
    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)
    
    du .= 0.0

    T = typeof(u.pressure[1])

    n_cells = length(cell_volumes)

    #grad_p = Vector{MVector{3, T}}(undef, n_cells)

    grad_p = [Real[0.0, 0.0, 0.0] for _ in 1:n_cells]
    
    # Loop over internal faces for Gradients
    for (i, (A, B)) in enumerate(cell_neighbor_map)
        area = connection_areas[i]
        norm = connection_normals[i] # Points A -> B
        
        p_face = 0.5 * (u.pressure[A] + u.pressure[B])
        
        contribution = p_face * area * norm
        #=
        if A == 1 && B == 2
            #println("grad p", grad_p[1])
            println(contribution)
        end
        =#
        grad_p[A] += contribution #this one, (2nd issue)
        grad_p[B] -= contribution #this one, (3rd issue)
    end

    # Loop over Boundary faces for Gradients (Crucial for Neumann/Dirichlet consistency)
    for (i, (cell_id, face_idx)) in enumerate(unconnected_cell_face_map)
        area = unconnected_areas[cell_id][face_idx] 
        norm = unconnected_normals[cell_id][face_idx]
        # Note: You need the normal for the unconnected face. 
        # Your current rebuild_fvm_geometry returns areas as scalars for unconnected?
        # You effectively need the normal here. 
        # For now, assuming Zero Gradient at boundaries for P reconstruction (Simplified):
        p_face = u.pressure[cell_id] 
        # We need the normal vector calculation for unconnected faces in your geometry struct.
        # skipping explicit boundary gradient addition for this snippet to keep it running, 
        # but technically required.
        contribution = p_face * area * norm
        grad_p[cell_id] += contribution #this one, (1st issue)
    end

    # Normalize Gradients by Volume
    for i in 1:n_cells
        grad_p[i] = grad_p[i] / u.pressure[i] #and this one, (4th issue)
        #all cause errors for ForwardDiff.gradient erroring with 
        #=
            LoadError: MethodError: no method matching Float64(::ForwardDiff.Dual{ForwardDiff.Tag{typeof(CostFun),Float64},Float64,3})
            Closest candidates are:
            Float64(::Real, ::RoundingMode) where T<:AbstractFloat at rounding.jl:200
            Float64(::T) where T<:Number at boot.jl:715
            Float64(::Int8) at float.jl:60
        =#
    end

    for (i, (A, B)) in enumerate(cell_neighbor_map)
        #A is the cell_id of the current cell and B is the cell_id of the neighboring cell
        mat_a = cell_mat_id_map[A]
        mat_b = cell_mat_id_map[B]

        area = connection_areas[i]
        norm = connection_normals[i]
        dist = connection_distances[i]

        rho_a = fluid_phys[mat_a].rho
        rho_b = fluid_phys[mat_b].rho
        rho_avg = 0.5 * (rho_a + rho_b) 

        mu_a = fluid_phys[mat_a].mu
        mu_b = fluid_phys[mat_b].mu
        mu_avg = 0.5 * (mu_a + mu_b) 

        avg_vel_x = 0.5 * (u.vel_x[A] + u.vel_x[B])
        avg_vel_y = 0.5 * (u.vel_y[A] + u.vel_y[B])
        avg_vel_z = 0.5 * (u.vel_z[A] + u.vel_z[B])

        vn_avg = (avg_vel_x * norm[1] + avg_vel_y * norm[2] + avg_vel_z * norm[3])

        p_diff = u.pressure[B] - u.pressure[A]

        vol_avg = 0.5 * (cell_volumes[A] + cell_volumes[B])
        rhie_chow_d = 0.5 * vol_avg

        grad_p_avg = 0.5 * (grad_p[A] + grad_p[B])
        grad_p_proj = dot(grad_p_avg, norm)

        grad_diff = (p_diff / dist) - grad_p_proj

        m_dot = rho_avg * area * (vn_avg - rhie_chow_d * grad_diff)

        du.pressure[A] -= m_dot
        du.pressure[B] += m_dot

        diff_flux_x = mu_avg * (u.vel_x[B] - u.vel_x[A]) / dist *  area
        diff_flux_y = mu_avg * (u.vel_y[B] - u.vel_y[A]) / dist *  area
        diff_flux_z = mu_avg * (u.vel_z[B] - u.vel_z[A]) / dist *  area

        #Convection momentum flux
        u_face_x = upwind(u.vel_x[A], u.vel_x[B], m_dot)
        u_face_y = upwind(u.vel_y[A], u.vel_y[B], m_dot)
        u_face_z = upwind(u.vel_z[A], u.vel_z[B], m_dot)
        
        conv_flux_x = m_dot * u_face_x
        conv_flux_y = m_dot * u_face_y
        conv_flux_z = m_dot * u_face_z

        #Pressure momentum flux
        p_face_val = 0.5 * (u.pressure[A] + u.pressure[B])
        press_x = p_face_val * area * norm[1]
        press_y = p_face_val * area * norm[2]
        press_z = p_face_val * area * norm[3]

        du.vel_x[A] += (diff_flux_x - conv_flux_x - press_x)
        du.vel_x[B] -= (diff_flux_x - conv_flux_x - press_x)
        
        du.vel_y[A] += (diff_flux_y - conv_flux_y - press_y)
        du.vel_y[B] -= (diff_flux_y - conv_flux_y - press_y)
        
        du.vel_z[A] += (diff_flux_z - conv_flux_z - press_z)
        du.vel_z[B] -= (diff_flux_z - conv_flux_z - press_z)
    end

    # Source and Capacity Loop
    #=
    for cell_id in bc_sys.free_idxs
        something here when it's needed
    end
    =#

    #=
    for cell_id in bc_sys.dirichlet_idxs
        du.temp1[cell_id] = 0.0
        du.temp2[cell_id] = 0.0
    end
    =#
end

grid_dimensions = (10, 5, 5)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

cell_dist = (right[1] / grid_dimensions[1])
cell_half_dist = (right[1] / grid_dimensions[1]) / 2 
addcellset!(grid, "inlet", x -> x[1] <= left[1] + cell_dist)
addcellset!(grid, "internal cells", x -> x[1] >= left[1] + cell_half_dist && x[1] <= right[1] - cell_half_dist)
addcellset!(grid, "outlet", x -> x[1] >= right[1] - cell_dist)

air_physics = FluidPhysics(401.0, 1.225, 1.8e-5, 1.0) #note that if we created a new struct for copper_physics performance would die
#steel_physics = FluidPhysics(30.0, 8000.0, 0.01, 1.0)

fluid_phys_vec = [air_physics, air_physics, air_physics]#, steel_physics]

inlet_cell_set_idxs = Set(getcellset(grid, "inlet"))
internal_cell_set_idxs = Set(getcellset(grid, "internal cells"))
outlet_cell_set_idxs = Set(getcellset(grid, "outlet"))

struct CellSet
    mat_id::Int
    cell_set_idxs::Set{Int}
end

inlet_set = CellSet(1, inlet_cell_set_idxs)
internal_cells_set = CellSet(2, internal_cell_set_idxs)
outlet_set = CellSet(3, outlet_cell_set_idxs)

cell_sets = [inlet_set, internal_cells_set, outlet_set]

free_idxs = Int[]
dirichlet_idxs = Int[]

#use :Dirichlet to fix a value to the initial value
function my_bc_mapper(cell_id)
    if cell_id in inlet_cell_set_idxs
        vel_x_bc = VelBC(:Dirichlet, 0.0) 
        vel_y_bc = VelBC(:Dirichlet, 0.0) 
        vel_z_bc = VelBC(:Dirichlet, 0.0) 
        pressure_bc = PressureBC(:Dirichlet, 200000)
        push!(dirichlet_idxs, cell_id)
        return [vel_x_bc, vel_y_bc, vel_z_bc, pressure_bc] 
    elseif cell_id in internal_cell_set_idxs
        vel_x_bc = VelBC(:Neumann, 0.0) 
        vel_y_bc = VelBC(:Neumann, 0.0) 
        vel_z_bc = VelBC(:Neumann, 0.0) 
        pressure_bc = PressureBC(:Neumann, 100000)
        push!(free_idxs, cell_id)
        return [vel_x_bc, vel_y_bc, vel_z_bc, pressure_bc] 
    elseif cell_id in outlet_cell_set_idxs
        vel_x_bc = VelBC(:Neumann, 0.0) 
        vel_y_bc = VelBC(:Neumann, 0.0) 
        vel_z_bc = VelBC(:Neumann, 0.0) 
        pressure_bc = PressureBC(:Neumann, 100000)
        push!(free_idxs, cell_id)
        return [vel_x_bc, vel_y_bc, vel_z_bc, pressure_bc] 
    end
end

vel_x_bcs = VelBC[]
vel_y_bcs = VelBC[]
vel_z_bcs = VelBC[]
pressure_bcs = PressureBC[]

n_cells = length(grid.cells)

for cell_id in 1:n_cells
    bcs = my_bc_mapper(cell_id) #returns vector [BC1, BC2, BC3, BC4]
    push!(vel_x_bcs, bcs[1])
    push!(vel_y_bcs, bcs[2])
    push!(vel_z_bcs, bcs[3])
    push!(pressure_bcs, bcs[4])
end

boundary_map = MultiPhysicsBCs(vel_x_bcs, vel_y_bcs, vel_z_bcs, pressure_bcs)

bc_sys = BoundarySystem(boundary_map, free_idxs, dirichlet_idxs)

db_grid_dimensions = (5, 5, 5)
db_left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
db_right = Ferrite.Vec{3}((1.1, 1.1, 1.1)) #has to be larger than the original grid
db_grid = generate_grid(Hexahedron, db_grid_dimensions, db_left, db_right)

u_proto = ComponentArray(vel_x = zeros(n_cells), vel_y = zeros(n_cells), vel_z = zeros(n_cells), pressure = zeros(n_cells))

for cell_id in 1:n_cells
    u_proto.vel_x[cell_id] = bc_sys.boundary_map.vel_x[cell_id].initial
    u_proto.vel_y[cell_id] = bc_sys.boundary_map.vel_y[cell_id].initial
    u_proto.vel_z[cell_id] = bc_sys.boundary_map.vel_z[cell_id].initial
    u_proto.pressure[cell_id] = bc_sys.boundary_map.pressure[cell_id].initial
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

cell_neighbor_map, neighbor_map_respective_node_ids = get_neighbor_map(grid)

unconnected_cell_face_map, unconnected_map_respective_node_ids = get_unconnected_map(grid)

nodes_of_cells = get_nodes_of_cells(grid)

proto_p_guess = get_node_coordinates(db_grid) #27 element Vector{Vector{}

p_guess = reduce(vcat, proto_p_guess) #81 element Vector{} # we do this because Optimization.solve(...) complains about vectors of vectors

db_nodes_of_cells = get_nodes_of_cells(db_grid)

cell_map, intrinsic_coordinates = precompute_ffd_data(proto_p_guess, db_nodes_of_cells, initial_node_coordinates)

db_grid_n_nodes = length(db_grid.nodes)

unflattened_p = eachcol(reshape(p_guess, 3, db_grid_n_nodes))

new_node_coordinates = apply_ffd_motion(unflattened_p, db_nodes_of_cells, cell_map, intrinsic_coordinates)


#= 
NOTE: We might have to declare these as constant to prevent dynamic dispatch in the future but we're leaving this out for now
const cell_neighbor_map = cell_neighbor_map 
const const_neighbor_map_respective_node_ids = neighbor_map_respective_node_ids
=#

f_closure = (du, u, p) -> FVM_iter_f!(
    du, u, p, 
    cell_neighbor_map, neighbor_map_respective_node_ids, 
    unconnected_cell_face_map, unconnected_map_respective_node_ids,
    nodes_of_cells, db_nodes_of_cells,
    cell_map, intrinsic_coordinates, 
    cell_mat_id_map, bc_sys, fluid_phys_vec, u_axes, db_grid_n_nodes
)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

du0 = u0 .* 0.0

f_closure(du0, u0, p_guess)

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
    cell_mat_id_map, bc_sys, fluid_phys_vec, u_axes, db_grid_n_nodes
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
    tmp_prob = remake(prob, u0=convert.(eltype(θ),prob.u0), p=θ) #not even sure this is required vs just passing in prob to sol = (...)
    sol = solve(tmp_prob, NonlinearSolve.NewtonRaphson(), p=θ, sensealg=SteadyStateAdjoint(diff_type = AutoEnzyme()), verbose = true, alias=SciMLBase.NonlinearAliasSpecifier(alias_u0=true))
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

    #current_vol = get_total_volume(θ) 
    #println(current_vol)
    #vol_loss = abs2(current_vol - initial_volume)

    #println("scaled_physics_loss, ", scaled_physics_loss)

    return (1.0 * scaled_physics_loss) #+ (vol_loss)
end

using ForwardDiff
ForwardDiff.gradient(loss, p_guess)
#This takes soooooooooo long

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
    callback = cb,
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

cell_volumes, 
cell_centroids, 
connection_areas, 
connection_normals, 
connection_distances, 
unconnected_areas = rebuild_fvm_geometry(cell_neighbor_map, neighbor_map_respective_node_ids, unconnected_cell_face_map, unconnected_map_respective_node_ids, new_node_coordinates, nodes_of_cells)


#=
TODO list
    - Use ComponentArrays to handle instances of metal wall/fluid cell interfaces as they are techically connected, but shouldn't be considered connected
        - We should probably make this a static topology like cell_neighbor_map that's passed into the function as constant once 
    - Implement FEM to calculate wall strain for the microchannel SMR
    - Deforming the mesh based on the FEM solve is not worth implementing at least for right now
=#

#=
What I have learned after a lot of debugging 
    - We should probably just give up on fixing zygote and ReverseDiff
    - Instead, we should focus on getting Enzyme to work by preventing any dynamic dispatch because it supports array mutation
    - While we try to get that to work, we should just use AutoForwardDiff as the deform box 

=#

#= Warning the code after this usually takes forever and errors 
function loss_enzyme(θ)
    pred_u = predict(θ)
    scaled_physics_loss = physics_loss(pred_u) / initial_physics_loss 
    #current_vol = get_total_volume(θ)
    #vol_loss = abs2(current_vol - initial_volume)
    return (1.0 * scaled_physics_loss) #+ vol_loss
end

# Use Enzyme directly instead of Zygote
d_p = Enzyme.make_zero(p_guess)
Enzyme.autodiff(Reverse, loss_enzyme, Active, Duplicated(p_guess, d_p))
try
    Enzyme.autodiff(Reverse, loss_enzyme, Active, Duplicated(p_guess, d_p))
    println("It Worked!")
catch err
    code_typed(err)
end
=#
#Tried doing this with an explicitly solved problem and I got the same error message, code_typed doesn't even do anything 

