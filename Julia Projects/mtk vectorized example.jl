using ModelingToolkit
using DifferentialEquations
using Ferrite
using SparseArrays
using LinearAlgebra

import ModelingToolkit: t_nounits as t, D_nounits as D

# --- Standard MTK Connectors ---
@connector HeatPort begin
    @variables begin
        T(t)
        Q_flow(t), [connect = Flow]
    end
end

function HeatConductionDomain(; name, grid, k_thermal, rho, cp)
    n_cells = getncells(grid)
    
    cell_centers = [sum(getcoordinates(cell))/length(getcoordinates(cell)) for cell in CellIterator(grid)]
    volumes = []
    cell_qr = QuadratureRule{RefHexahedron}(2)
    cell_values = CellValues(cell_qr, Lagrange{RefHexahedron, 1}())
    
    for cell in CellIterator(grid)
        vol = 0.0 
        for qp in 1:getnquadpoints(cell_values)
            d_vol = getdetJdV(cell_values, qp) 
            vol += d_vol
        end
        push!(volumes, vol)
    end
    
    # Build Sparse Conductance Matrix (K)
    # Row i, Col j = Conductance between i and j.
    # Row i, Col i = Sum of conductances (negative)
    I = Int[]
    J = Int[]
    V = Float64[]
    
    topology = ExclusiveTopology(grid)
    
    for i in 1:n_cells
        neighbors = [] 
        for face_idx in 1:nfacets(grid.cells[i])
             neighbor_info = topology.face_face_neighbor[i, face_idx]
             if !isempty(neighbor_info)
                 n_idx = collect(neighbor_info[1].idx)[1] # Neighbor ID
                 
                 # Calculate geometric conductance
                 dist = norm(cell_centers[i] - cell_centers[n_idx])
                 # Approximate face area (hack for hex)
                 area = (volumes[i])^(2/3) 
                 G = k_thermal * area / dist
                 
                 push!(I, i); push!(J, n_idx); push!(V, G)
                 push!(I, i); push!(J, i);     push!(V, -G)
             end
        end
    end
    
    K_mat = sparse(I, J, V, n_cells, n_cells)
    C_vec = volumes .* rho .* cp
    
    @variables T(t)[1:n_cells] = [300.0 for i in 1:n_cells]
    @parameters C[1:n_cells] = C_vec
    
    left_set = getcellset(grid, "left")
    right_set = getcellset(grid, "right")
    
    left_idxs = collect(left_set)
    right_idxs = collect(right_set)
    
    port_left = HeatPort(name=:port_left)
    port_right = HeatPort(name=:port_right)
    
    # Vectorized PDE: C * dT/dt = K * T + Q_external
    
    eqs = Equation[]
    
    conduction_term = K_mat * T 
    
    for i in 1:n_cells
        rhs = conduction_term[i]
        
        if i in left_idxs
            rhs += port_left.Q_flow / length(left_idxs)
        end
        
        if i in right_idxs
            rhs += port_right.Q_flow / length(right_idxs)
        end
        
        push!(eqs, D(T[i]) ~ rhs / C[i])
    end

    #average T on boundary face for ports of bcs components
    
    T_left = 0
    for left_idx in left_idxs 
        T_left += T[left_idx]
    end
    T_right = 0
    for right_idx in right_idxs 
        T_right += T[right_idx]
    end

    push!(eqs, port_left.T ~ T_left / length(left_idxs))
    push!(eqs, port_right.T ~ T_right / length(right_idxs))

    return ODESystem(eqs, t, [T...], [C...]; systems=[port_left, port_right], name=name)
end

#generate mesh
left_coord = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right_coord = Ferrite.Vec{3}((1.0, 1.0, 1.0))

grid_dim = (2, 2, 2) 

grid = generate_grid(Hexahedron, grid_dim, left_coord, right_coord)


length_to_node_ratio = right_coord[1] / collect(grid_dim)[1]

addcellset!(grid, "left", x -> x[1] <= left_coord[1] + length_to_node_ratio)
getcellset(grid, "left")
addcellset!(grid, "right", (x) -> x[1] >= (right_coord[1] - 0.0000001) - (length_to_node_ratio)) #1 doesn't work
getcellset(grid, "right")

@named pipe = HeatConductionDomain(grid=grid, k_thermal=200, rho=2700, cp=900)

@mtkmodel ConstantTemperature begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        T_fixed = 350.0
    end
    @equations begin
        port.T ~ T_fixed
    end
end

@mtkmodel ConstantHeatFlow begin
    @components begin
        port = HeatPort()
    end
    @parameters begin
        Q_fixed = 5000.0
    end
    @equations begin
        port.Q_flow ~ -Q_fixed # Negative because flow leaves the source
    end
end

@named inlet_bc = ConstantHeatFlow(Q_fixed=10000.0)
#@named outlet_bc = ConstantHeatFlow(Q_fixed=10000.0) #ConstantTemperature Doesn't work here for some reason
#ConstantHeatFlow(Q_fixed=10000.0)
#ConstantTemperature(T_fixed=300.0)

connections = [
    connect(inlet_bc.port, pipe.port_left)
    #connect(pipe.port_right, outlet_bc.port)
]

@mtkcompile model = ODESystem(connections, t, systems=[pipe, inlet_bc, #=outlet_bc=#]) #mtk compile seems significantly faster than structural_simplify

equations(expand_connections(model))

#simplify
#@time sys = structural_simplify(model) 

prob = ODEProblem(model, [], (0.0, 10.0))
@time sol = solve(prob, FBDF()) 

#using Plots
#center_idx = Int(round(getncells(grid)/2))
#plot(sol, idxs=[pipe.T[center_idx], pipe.port_left.T])