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