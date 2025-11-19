#Check out FiniteVolumeMethod.jl, specifically this function: FVMGeometry(gmsh.model)
#We should probably just use FiniteVolumeMethod solvers directly instead of trying to couple it with mtkcomponents
using Ferrite
using SparseArrays

grid = generate_grid(Quadrilateral, (2, 2))

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
#START of unique vertex connections finder for FEM
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
#END of unique vertex connections finder for FEM


#oh, VertexIndex(x, y) is where x is the cell and y is the vertex


#START of FVM attempt
#START of face area finder
left = Vec{3}((-1.0, -1.0, -1.0))
right = Vec{3}((1.0, 1.0, 1.0))

grid = generate_grid(Hexahedron, (2, 2, 2), left, right)

qr = FacetQuadratureRule{RefHexahedron}(2) #no idea what the 2 represents here #represents the number of integration points, basically higher number = higher accuracy but more computation

ip = Lagrange{RefHexahedron, 1}() #no idea what the 1 represents here #1 = linear elements 2 = quadratic/curved edges

fv = FacetValues(qr, ip)

cell_coords = getcoordinates(grid, 1)

reinit!(fv, cell_coords, 4)

getdetJdV(fv, 2) #not sure what the 2 represents here #The 2 is the quadrature point index.

area = []
for i in 1:getnquadpoints(fv)
    push!(area, getdetJdV(fv, i))
end
area
#END of face area finder

#START of cell volume finder
ip = Lagrange{RefHexahedron, 1}() 
qr = QuadratureRule{RefHexahedron}(2) 
cellvalues = CellValues(qr, ip) #we get CellValues instead of FacetValues here

vol_per_cell = []

for cc in CellIterator(grid)
    reinit!(cellvalues, cc)
    vol_i = 0.0 # Initialize cell volume/area
    for qp in 1:getnquadpoints(cellvalues)
        d_vol = getdetJdV(cellvalues, qp) # Get the integration weight (which includes the Jacobian and quadrature weight)
        vol_i += d_vol
    end
    push!(vol_per_cell, vol_i)
end
return vol_per_cell
#END of cell volume finder

#START of combined loop
left = Vec{3}((0.0, 0.0, 0.0))
right = Vec{3}((1.0, 1.0, 1.0))

grid = generate_grid(Hexahedron, (2, 2, 2), left, right)

ip = Lagrange{RefHexahedron, 1}() #no idea what the 1 represents here

cell_qr = QuadratureRule{RefHexahedron}(2)
cell_values = CellValues(cell_qr, ip)

facet_qr = FacetQuadratureRule{RefHexahedron}(2) #no idea what the 2 represents here
facet_values = FacetValues(facet_qr, ip)

cell_coords = getcoordinates(grid, 1) #this returns a cube with 1.0 long sides 
#reinit!(facet_values, cell_coords, 1) #you get NaN for neighbor area if you don't do this

Ferrite.facets(getcells(grid)[1]) 

cell_dimensions = []

for cell in CellIterator(grid)
    reinit!(cell_values, cell)
    vol_i = 0.0 
    for qp in 1:getnquadpoints(cell_values)
        d_vol = getdetJdV(cell_values, qp)
        vol_i += d_vol
    end
    
    area = []
    cell_coords = getcoordinates(cell)
    
    for facet_index in 1:nfacets(cell)  
        reinit!(facet_values, cell_coords, facet_index)
        
        area_i = 0.0
        for qp in 1:getnquadpoints(facet_values)
            area_i += getdetJdV(facet_values, qp)
        end
        
        push!(area, area_i)
    end
    push!(cell_dimensions, [[area], vol_i])
end
cell_dimensions 

#END of combined loop
#END of FVM attempt
#Don't know if these would be useful 
struct VertexConnection
    vertex1_idx::Int
    vertex2_idx::Int
    coords1::Vec{2, Float64}  #or Vec{3, Float64} for 3D
    coords2::Vec{2, Float64}
    distance::Float64
    cross_section_area::Float64  #For 2D: element thickness × perpendicular width
end

struct FacetConnection
    cell1_idx::Int
    cell2_idx::Int
    facet_area::Float64
    facet_normal::Vec{2, Float64}  # Points from cell1 to cell2
    distance_between_centroids::Float64
end

#mtk component matching
struct FEMNodeMapping{T}
    vertex_idx::Int
    component::T  #ex. HeatCapacitor
end

struct FEMEdgeMapping{T}
    connection::VertexConnection
    component::T  #ex. ThermalConductor
end

struct FVMCellMapping{T}
    cell_idx::Int
    component::T
    volume::Float64
end

struct FVMFacetMapping{T}
    connection::FacetConnection
    component::T  #ex. some flux component
end



"""

#oh, VertexIndex(x, y) is where x is the cell and y is the vertex
qr = QuadratureRule{RefTriangle}(1)

FacetQuadratureRule()

ip = Lagrange{RefTriangle, 1}()

FacetValues(qr, ip)
#Creates an object for calculations on cell facets (edges in 2D, faces in 3D)
#Stores quadrature rules and shape functions for facet integration

FaceVectorValues() 
#Same as FacetValues but for vector-valued fields
#Use when interpolating vector quantities (displacement, velocity, etc.)

getdetJdV()
#Summing getdetJdV over all qp gives the area/volume

getcoordinates() #returns the physical coordinates of every vertex in a cell
"""




