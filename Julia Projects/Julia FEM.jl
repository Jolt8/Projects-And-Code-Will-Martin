#Check out FiniteVolumeMethod.jl, specifically this function: FVMGeometry(gmsh.model)
#We should probably just use FiniteVolumeMethod solvers directly instead of trying to couple it with mtkcomponents
using Ferrite
using SparseArrays

grid = generate_grid(Quadrilateral, (1, 1))

getcells(grid)

getnodes(grid)

Ferrite.nnodes_per_cell(grid, 1) #you can get the amount of nodes per cell for a specific index by using this

cell_pos = 4

Ferrite.vertices(getcells(grid)[cell_pos])

Ferrite.edges(getcells(grid)[cell_pos])

Ferrite.faces(getcells(grid)[cell_pos])

Ferrite.facets(getcells(grid)[cell_pos])

node_pos = 2

get_node_coordinate(getnodes(grid)[node_pos])

ExclusiveTopology

top = ExclusiveTopology(grid)

top.face_face_neighbor

top.vertex_vertex_neighbor

top.edge_edge_neighbor

Ferrite.vertex_star_stencils(top, grid)

Ferrite.getspatialdim(grid) #ooh, this is useful, it tells you whether or not the grid is 1D, 2D, or 3D

VertexIndex

Ferrite.getneighborhood(top, grid, VertexIndex((1, 1)))

Ferrite.getneighborhood(top, grid, EdgeIndex((1, 2)))

#below is what we need
top = ExclusiveTopology(grid)

vertex_neighbors = Ferrite.vertex_star_stencils(top, grid)

unique_connections = Set{Tuple{Int, Int}}()

for i in eachindex(vertex_neighbors)
    for j in 2 #1, 1 #starts at 2: to not include (1, 1)
        println(i, ", ", collect(vertex_neighbors[i][j].idx)[2])
        push!(unique_connections, (min(i, collect(vertex_neighbors[i][j].idx)[2]), max(i, collect(vertex_neighbors[i][j].idx)[2])))
    end
end

unique_connections = collect.(unique_connections)


#oh, VertexIndex(x, y) is where x is the cell and y is the vertex


