using Ferrite
using SparseArrays

function grid_vertices_unique_connections(grid)
    top = ExclusiveTopology(grid)

    vertex_neighbors = Ferrite.vertex_star_stencils(top, grid)

    unique_connections = Set{Tuple{Int, Int}}()

    for i in eachindex(vertex_neighbors)
        node_at_i = collect(vertex_neighbors[i][1].idx)[2]
        for j in 2:length(vertex_neighbors[i])
            neighbor_at_i = collect(vertex_neighbors[i][j].idx)[2]
            push!(unique_connections, (min(node_at_i, neighbor_at_i), max(node_at_i, neighbor_at_i)))
        end
    end
    return collect(unique_connections)
end

left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid = generate_grid(Hexahedron, (1, 1, 1), left, right)

grid_vertices_unique_connections(grid) #1 cell with 8 vertices gives 12 unique connections so it works!

grid_2D = generate_grid(Quadrilateral, (1, 1))

grid_vertices_unique_connections(grid_2D) #gives 4 unique connections so it also works!
