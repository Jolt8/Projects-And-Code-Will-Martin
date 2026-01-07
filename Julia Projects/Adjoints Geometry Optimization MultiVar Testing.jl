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



# --- 1. Data Structures ---
# We strip physics values out of the static structs

struct CellData
    volume::Float64
    centroid::Vec{3, Float64}
    material_id::Int # Added this to map parameters to cells
end

struct Connection
    cell_idx_a::Int
    cell_idx_b::Int
    area::Float64
    normal::Vec{3, Float64}
    distance::Float64
end

struct HeatBC
    type::Symbol 
    initial::Vector{Float64}
    # Physics is removed from here, determined by material_id at runtime
end

struct FVMGeometry
    mesh_cells::Vector{CellData}
    connections::Vector{Connection}
    boundary_map::Vector{HeatBC}
    free_idxs::Vector{Int} 
    dirichlet_idxs::Vector{Int}
end

struct PhysicsParams
    k_cu::Float64; rho_cu::Float64; cp_cu::Float64
    k_st::Float64; rho_st::Float64; cp_st::Float64
end

# --- 2. Physics Helper Functions ---
# Modified to take raw numbers (k, rho, cp) instead of a struct

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

# --- 3. Grid Generation & Setup ---

grid_dimensions = (5, 3, 3) # Smaller grid for testing AD speed
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

cell_half_dist = (right[1] / grid_dimensions[1]) / 2 
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))

copper_cell_set_idxs = Set(getcellset(grid, "copper"))
steel_cell_set_idxs = Set(getcellset(grid, "steel"))

cell_sets = [copper_cell_set_idxs, steel_cell_set_idxs]

function build_fvm_geometry(grid, cell_sets)
    n_cells = getncells(grid)
    ref_shape = getrefshape(grid.cells[1])
    poly_interp = Lagrange{ref_shape, 1}()
    facet_qr = FacetQuadratureRule{ref_shape}(2)
    facet_values = FacetValues(facet_qr, poly_interp)
    cell_qr = QuadratureRule{ref_shape}(2)
    cell_values = CellValues(cell_qr, poly_interp)

    cells_data = Vector{CellData}(undef, n_cells)
    
    # 1. Build Cell Data with Material IDs
    for cell in CellIterator(grid)
        id = cellid(cell)
        Ferrite.reinit!(cell_values, cell)
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values))
        coords = getcoordinates(cell)
        cent = sum(coords) / length(coords)
        
        # Assign Material ID (1 for Copper, 2 for Steel)
        mat_id = (id in cell_sets[1]) ? 1 : 2
        cells_data[id] = CellData(vol, cent, mat_id)
    end

    top = ExclusiveTopology(grid)
    connections = Vector{Connection}()
    boundary_map = Vector{HeatBC}(undef, n_cells)
    free_idxs = Vector{Int}()
    dirichlet_idxs = Vector{Int}()
    
    # 2. Build Connections and BCs
    for i in 1:n_cells
        # Simplified BC logic for demonstration
        if i in cell_sets[1]
            boundary_map[i] = HeatBC(:Free, [1000.0, 700.0])
        else
            boundary_map[i] = HeatBC(:Free, [500.0, 300.0])
        end

        if boundary_map[i].type == :Dirichlet
            push!(dirichlet_idxs, i)
        else
            push!(free_idxs, i)
        end

        for face_idx in 1:nfacets(grid.cells[i])
            neighbor_info = top.face_face_neighbor[i, face_idx]
            if !isempty(neighbor_info)
                neighbor_idx = collect(neighbor_info[1].idx)[1]
                if i < neighbor_idx
                    coords = getcoordinates(grid, i)
                    Ferrite.reinit!(facet_values, coords, face_idx)
                    area = sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values))
                    n_ref = getnormal(facet_values, 1) 
                    centroid_a = cells_data[i].centroid
                    centroid_b = cells_data[neighbor_idx].centroid
                    dist = norm(centroid_b - centroid_a)
                    push!(connections, Connection(i, neighbor_idx, area, n_ref, dist))
                end
            end
        end
    end
    
    return FVMGeometry(cells_data, connections, boundary_map, free_idxs, dirichlet_idxs)
end

println("Building FVM Geometry...")
fvm_geo = build_fvm_geometry(grid, cell_sets)

#=
struct CellData
    volume::Float64
    centroid::Vec{3, Float64}
    material_id::Int # Added this to map parameters to cells
end

struct Connection
    cell_idx_a::Int
    cell_idx_b::Int
    area::Float64
    normal::Vec{3, Float64}
    distance::Float64
end
=#


fvm_geo.connections

grid_dimensions

struct GridTopology
    # Map: Cell Index -> List of 8 Node Indices (for Hexahedrons)
    cell_node_indices::Vector{Vector{Int}} 
    
    # Map: Connection Index -> List of 4 Node Indices (for the shared face)
    conn_face_indices::Vector{Vector{Int}}
    
    # Base coordinates (undeformed) to apply deltas to
    base_coords::Vector{Vec{3, Float64}} 
    
    # Map: Cell Index -> (i, j, k) indices (useful for applying p_shape)
    cell_ijk::Vector{Tuple{Int, Int, Int}}
    
    # Map: Node Index -> (i, j, k) indices (useful for deforming nodes)
    node_ijk::Vector{Tuple{Int, Int, Int}}
end 

grid.nodes[3]

getcoordinates(grid.nodes[3])

get_node_coordinate(grid.nodes[3])

function build_grid_topology(grid, fvm_geo)
    n_cells = getncells(grid)
    n_nodes = getnnodes(grid)
    
    # 1. Extract Base Coords
    base_coords = [get_node_coordinate(grid.nodes[i]) for i in 1:n_nodes]

    # 2. Map Cells to Nodes
    cell_node_indices = Vector{Vector{Int}}(undef, n_cells)
    for i in 1:n_cells
        cell_node_indices[i] = [grid.cells[i].nodes...] # Unpack tuple to vector
    end
    
    # 3. Map Connections to Face Nodes
    # This is the tricky part: finding which nodes are shared
    conn_face_indices = Vector{Vector{Int}}(undef, length(fvm_geo.connections))
    
    for (k, conn) in enumerate(fvm_geo.connections)
        nodes_a = Set(cell_node_indices[conn.cell_idx_a])
        nodes_b = Set(cell_node_indices[conn.cell_idx_b])
        # The shared face is the intersection of nodes
        shared_nodes = collect(intersect(nodes_a, nodes_b))
        
        # For a Hex grid, a face MUST have 4 nodes
        @assert length(shared_nodes) == 4 "Connection $k has $(length(shared_nodes)) shared nodes, expected 4."
        conn_face_indices[k] = shared_nodes
    end

    # 4. Generate (i,j,k) lookups (Assuming structured grid generated by Ferrite)
    # Note: This relies on Ferrite's specific node ordering for structured grids.
    # For a grid_dimensions = (Nx, Ny, Nz)
    Nx, Ny, Nz = grid_dimensions
    
    # Ferrite generation usually orders x fast, then y, then z.
    # We can reverse engineer i,j,k from the centroid or just trust the ordering if standard.
    # Let's use the centroids to be robust.
    
    node_ijk = Vector{Tuple{Int,Int,Int}}(undef, n_nodes)
    # Simple bounding box logic to bucket nodes
    xs = sort(unique([c[1] for c in base_coords]))
    ys = sort(unique([c[2] for c in base_coords]))
    zs = sort(unique([c[3] for c in base_coords]))
    
    for n in 1:n_nodes
        c = base_coords[n]
        # Find index in sorted list (inefficient but run-once)
        i = findfirst(x->isapprox(x, c[1], atol=1e-5), xs)
        j = findfirst(x->isapprox(x, c[2], atol=1e-5), ys)
        k = findfirst(x->isapprox(x, c[3], atol=1e-5), zs)
        node_ijk[n] = (i, j, k)
    end
    
    # Do same for cells (optional, but useful for mapping p_shape to cells)
    cell_ijk = Vector{Tuple{Int,Int,Int}}(undef, n_cells)
    # ... (similar logic using fvm_geo.mesh_cells[i].centroid) ...
    # Placeholder for brevity:
    cell_ijk = [(1,1,1) for _ in 1:n_cells] 
    
    return GridTopology(cell_node_indices, conn_face_indices, base_coords, cell_ijk, node_ijk)
end

function hex_volume(nodes::Vector{<:AbstractVector})
    # Break hex into 5 tetrahedrons or use mean centroid method.
    # Robust method: Sum of 6 pyramids formed by faces and the centroid.
    
    # Calculate geometric centroid
    c = sum(nodes) / 8.0
    
    vol = 0.0
    # Face 1: 1-2-6-5 (Bottom)
    vol += pyr_vol(c, nodes[1], nodes[2], nodes[6], nodes[5])
    # Face 2: 2-3-7-6 (Right)
    vol += pyr_vol(c, nodes[2], nodes[3], nodes[7], nodes[6])
    # Face 3: 3-4-8-7 (Top)
    vol += pyr_vol(c, nodes[3], nodes[4], nodes[8], nodes[7])
    # Face 4: 4-1-5-8 (Left)
    vol += pyr_vol(c, nodes[4], nodes[1], nodes[5], nodes[8])
    # Face 5: 1-2-3-4 (Back) - Check winding order based on Ferrite standard
    vol += pyr_vol(c, nodes[1], nodes[4], nodes[3], nodes[2])
    # Face 6: 5-6-7-8 (Front)
    vol += pyr_vol(c, nodes[5], nodes[6], nodes[7], nodes[8])
    
    return vol
end

function pyr_vol(apex, p1, p2, p3, p4)
    # Volume of pyramid with square base = sum of two tetrahedrons
    # Tet volume = |(a-d) . ((b-d) x (c-d))| / 6
    v1 = det([p1-apex p2-apex p3-apex]) / 6.0
    v2 = det([p1-apex p3-apex p4-apex]) / 6.0
    return abs(v1) + abs(v2)
end

function face_area(p1, p2, p3, p4)
    # Area of a quad (split into 2 triangles)
    # tri1: 1-2-3, tri2: 1-3-4
    a1 = 0.5 * norm(cross(p2-p1, p3-p1))
    a2 = 0.5 * norm(cross(p3-p1, p4-p1))
    return a1 + a2
end

# 3b. The Geometry Update Function
# This function is fully differentiable.
# p_shape: Vector of floats controlling width at each X-index

function update_geometry(p_shape, topo::GridTopology)
    n_nodes = length(topo.base_coords)
    n_cells = length(topo.cell_node_indices)
    n_conns = length(topo.conn_face_indices)
    
    # --- 1. Deform Nodes ---
    # Create a new array of coordinates. 
    # Important: Do NOT mutate topo.base_coords. Create new array.
    
    # We map p_shape indices to the grid's I-index (x-axis)
    # p_shape[i] determines the multiplier for the Y-coordinate (Width)
    
    current_coords = Vector{Vec{3, eltype(p_shape)}}(undef, n_nodes)
    
    for n in 1:n_nodes
        idx_i, idx_j, idx_k = topo.node_ijk[n]
        
        # Get deformation factor for this X-slice
        # Clamp index just in case
        safe_i = min(max(idx_i, 1), length(p_shape))
        factor = p_shape[safe_i] 
        
        base = topo.base_coords[n]
        
        # Example Deformation: Vary Y (width) based on X position
        # factor = 1.0 means original width
        # factor = 2.0 means double width
        new_x = base[1]
        new_y = base[2] * factor 
        new_z = base[3]
        
        current_coords[n] = Vec{3}((new_x, new_y, new_z))
    end
    
    # --- 2. Calculate Cell Volumes & Centroids ---
    vols = Vector{eltype(p_shape)}(undef, n_cells)
    centroids = Vector{Vec{3, eltype(p_shape)}}(undef, n_cells)
    
    for i in 1:n_cells
        node_ids = topo.cell_node_indices[i]
        # Gather the 8 coordinates for this cell
        cell_nodes = [current_coords[nid] for nid in node_ids]
        
        vols[i] = hex_volume(cell_nodes)
        centroids[i] = sum(cell_nodes) / 8.0
    end
    
    # --- 3. Calculate Connection Areas & Distances ---
    areas = Vector{eltype(p_shape)}(undef, n_conns)
    dists = Vector{eltype(p_shape)}(undef, n_conns)
    
    for k in 1:n_conns
        # Distances
        idx_a = fvm_geo.connections[k].cell_idx_a # Access global constant or pass as arg? 
        # Ideally pass fvm_geo structure, but for now we assume indices line up
        idx_b = fvm_geo.connections[k].cell_idx_b
        
        c_a = centroids[idx_a]
        c_b = centroids[idx_b]
        dists[k] = norm(c_b - c_a)
        
        # Areas
        face_node_ids = topo.conn_face_indices[k]
        f_nodes = [current_coords[nid] for nid in face_node_ids]
        areas[k] = face_area(f_nodes[1], f_nodes[2], f_nodes[3], f_nodes[4])
    end
    
    return vols, centroids, areas, dists
end


topo = build_grid_topology(grid, fvm_geo)


p_shape = zeros(grid_dimensions[1]) .+= 1.0

p_shape[3] = 3.0

p_shape[5] = 5.0

update_geometry(p_shape, topo)


function FVM_iter_f!(du, u, p_params, 
                     geo_static::FVMGeometry, topo::GridTopology, 
                     phys::PhysicsParams)
    vols, centroids, areas, dists = update_geometry(p_params, topo)
    
    du .= 0.0
    
    # 2. Connection Loop
    for (k, conn) in enumerate(geo_static.connections)
        idx_a = conn.cell_idx_a
        idx_b = conn.cell_idx_b
        
        # Get Dynamic Geometric Data
        area = areas[k]
        dist = dists[k]
        
        # Physics ...
        k_avg = 400.0 # Example
        u_a = u[idx_a]
        u_b = u[idx_b]
        
        flux = numerical_flux(k_avg, u_a, u_b, area, dist)
        
        du[idx_a] -= flux
        du[idx_b] += flux
    end
    
    # 3. Source Loop
    for i in geo_static.free_idxs
        vol = vols[i] # Dynamic Volume
        # ... Capacity/Source logic using vol ...
        
    end
end
