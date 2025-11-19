left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid = generate_grid(Hexahedron, (2, 2, 2), left, right)

top = ExclusiveTopology(grid)

getncells(grid)

top.cell_neighbor

top.face_face_neighbor

top.face_face_neighbor[1, 5]
top.face_face_neighbor[1, 3]
top.face_face_neighbor[2, end] 
#If cell 1 had connection with cell 2 at (cell 1, face 3) and (cell 2, face 5)
#OH, I get how this is formatted, using the example above, 1 is the first index of the matrix for cell 1, 
#the 3rd index of the vector within the array is the index of the corresponding face idx for cell 1,
#that vector at the index [1, 3] of the matrix contains the struct FaceIndex that contains (2, 5)
#this means that cell 2 has a connection with cell 1 at cell 2's 5th side 

#this is confirmed by the fact that at index 5 (5th side) of the 2nd index (cell 2) in the matrix contains (FaceIndex(1, 3) or a connection with face 3 of cell 1

unique_connections = Set{Tuple{Tuple{Number, Number}, Tuple{Number, Number}}}() #we should use a set but I don't know how to set it up right

top.face_face_neighbor[1, 3]

for cell in CellIterator(grid)
    i = cellid(cell) 
    ff_neighbor_matrix = top.face_face_neighbor
    for j in 1:nfacets(cell)
        if isempty(ff_neighbor_matrix[i, j])
            nothing
        else
            #[i, j] = [1, 3] #[cell, face]
            #ff_neighbor_matrix[i, j] = [FacetIndex(2, 5)]
            #ff_neighbor_matrix[i, j][1].idx = (2, 5)
            #collect(ff_neighbor_matrix[i, j][1].idx) = [2, 5] #[cell, face]
            #neighbor_cell_and_face = collect(ff_neighbor_matrix[i, j][1].idx)
            neighbor_cell_idx = collect(ff_neighbor_matrix[i, j][1].idx)[1] #I don't know which of these is clearer (top vs bottom two)
            connecting_face_idx = collect(ff_neighbor_matrix[i, j][1].idx)[2]

            connection = minmax((i, j), (neighbor_cell_idx, connecting_face_idx))
            push!(unique_connections, connection)
        end
    end
end
unique_connections