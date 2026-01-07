using Ferrite
using DifferentialEquations
using LinearAlgebra
using Dates
using SparseArrays
using SciMLSensitivity
using Optimization, OptimizationPolyalgorithms, Zygote
using Enzyme

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
    initial::Float64
    # Physics is removed from here, determined by material_id at runtime
end

struct FVMGeometry
    mesh_cells::Vector{CellData}
    connections::Vector{Connection}
    boundary_map::Vector{HeatBC}
    free_idxs::Vector{Int} 
    dirichlet_idxs::Vector{Int}
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

grid_dimensions = (3, 3, 3) # Smaller grid for testing AD speed
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

cell_half_dist = (right[1] / grid_dimensions[1]) / 2 
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))

global copper_cell_set_idxs = Set(getcellset(grid, "copper"))
global steel_cell_set_idxs = Set(getcellset(grid, "steel"))

function build_fvm_geometry(grid)
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
        mat_id = (id in copper_cell_set_idxs) ? 1 : 2
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
        if i in copper_cell_set_idxs
            boundary_map[i] = HeatBC(:Neumann, 500.0)
        else
            boundary_map[i] = HeatBC(:Neumann, 300.0)
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
fvm_geo = build_fvm_geometry(grid)

# --- 4. The Parameter Vector ---
# p[1:3] = Copper (k, rho, cp)
# p[4:6] = Steel (k, rho, cp)
# We flatten this so Zygote can differentiate it easily.
p_true = [401.0, 8960.0, 385.0, 30.0, 8000.0, 460.0]
# Add source terms to p if you want to optimize them, otherwise hardcode or use aux data

# --- 5. The ODE Function (Closure) ---
# We pass fvm_geo as a captured variable (closure). 
# It is NOT in 'p', so AD ignores it (which is good).
function assemble_system(p, geo)
    # ... unpack p ...
    I, J, V = Int[], Int[], Float64[]
    
    for conn in geo.connections
        # Calculate conductance (transmissibility)
        k1, rho1, cp1 = p[1], p[2], p[3]
        k2, rho2, cp2 = p[4], p[5], p[6]

        mat_a = geo.mesh_cells[conn.cell_idx_a].material_id
        mat_b = geo.mesh_cells[conn.cell_idx_b].material_id

        # Get k based on material ID
        k_a = (mat_a == 1) ? k1 : k2
        k_b = (mat_b == 1) ? k1 : k2

        k_avg = 2 * k_a * k_b / (k_a + k_b)
        
        trans = k_avg * conn.area / conn.distance
        
        # Diagonal elements (loss)
        push!(I, conn.cell_idx_a); push!(J, conn.cell_idx_a); push!(V, -trans)
        push!(I, conn.cell_idx_b); push!(J, conn.cell_idx_b); push!(V, -trans)
        
        # Off-diagonal elements (gain from neighbor)
        push!(I, conn.cell_idx_a); push!(J, conn.cell_idx_b); push!(V, trans)
        push!(I, conn.cell_idx_b); push!(J, conn.cell_idx_a); push!(V, trans)
    end
    
    return sparse(I, J, V)
end

function FVM_Matrix_f!(du, u, p, t, geo)
    # This might be re-calculated every step if K depends on T (nonlinear)
    # Or just once if linear
    K = assemble_system(p, geo) 
    
    mul!(du, K, u) 
end

# Create the ODEProblem using a closure to capture `fvm_geo`
# The syntax is f(du, u, p, t)
f_closure = (du, u, p, t) -> FVM_Matrix_f!(du, u, p, t, fvm_geo)


u0 = assemble_system(p_true, fvm_geo)
u0 = zeros(n_cells, n_cells)
for i in 1:n_cells
    u0[i, :] .= fvm_geo.boundary_map[i].initial
end

u0


t0, tMax = 0.0, 1000.0 
desired_steps = 100
dt = tMax / 100
tspan = (t0, tMax)
t = t0:dt:tMax;
prob = ODEProblem(f_closure, u0, tspan, p_true)

desired_amount_of_u = 100

println("Solving Forward Problem...")
# Using Tsit5 for explicit, ensure dt is small enough for stability
@time sol = solve(prob, FBDF(), saveat=t)
stats = @timed sol = solve(prob, FBDF(), saveat=t)
solve_time = stats.time
@time sol = solve(prob, FBDF(), saveat=t)
arr_sol = Array(sol)[1, :, :] # The "Truth" data to train against


# Start with wrong parameters

p_guess = [350.0, 50.0] 
p_fixed = [8960.0, 385.0, 8000.0, 460.0]

function assemble_p(θ, constants)
    k1, k2 = θ[1], θ[2]
    rho1, cp1, rho2, cp2 = constants[1], constants[2], constants[3], constants[4]
    
    # Reconstruct the order expected by FVM_iter_f!
    # [k1, rho1, cp1, k2, rho2, cp2]
    return [k1, rho1, cp1, k2, rho2, cp2]
end

function predict(θ, constants)
    # solve needs to know to differentiate through θ (which enters as p)
    # sensealg=InterpolatingAdjoint() is standard for memory efficiency
    p_full = assemble_p(θ, constants)
    Array(solve(prob, FBDF(), p=p_full, saveat=t, sensealg=InterpolatingAdjoint(autodiff=AutoEnzyme())))[1, :, :]
end

zeros(Float64, length(arr_sol[:, 1])) .+ 400.0
#=
function loss(θ)
    pred = predict(θ)
    # Simple MSE
    return sum(abs2, pred .- (zeros(Float64, length(arr_sol[:, 1])) .+ 400.0))
end=#

test_pred = predict(p_guess, p_fixed)[4, :]

sum(abs2, test_pred[4, :] .- 400.0)

length(test_pred)

(sum(abs2, test_pred[4, :] .- 400.0) / length(test_pred))

function loss(θ, constants)
    pred = predict(θ, constants)
    # Simple MSE
    #return (sum(abs2, pred[4, :] .- 400.0) / length(test_pred))
    return sum(abs2, pred[4, :] .- 400.0) 
end

l = loss(p_guess, p_fixed)

LOSS = [] # Loss accumulator
PRED = [] # prediction accumulator
PARS = [] # parameters accumulator

cb = function (state, l) #callback function to observe training
    display(l)
    display(state.u)
    #pred = predict(state.u)
    #append!(PRED, [pred])
    append!(LOSS, l)
    append!(PARS, [state.u])
    false
end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction(loss, adtype)
optprob = Optimization.OptimizationProblem(optf, p_guess, p_fixed)

maxiters = 20

solve_time

max_time_to_solve = maxiters * solve_time #this is usually wrong, it's often much greater 

println("Starting Optimization...")

proceed = true

if proceed 
    res = Optimization.solve(optprob, PolyOpt(), callback=cb, maxiters=maxiters)
end

PARS[end]

#PARS[end] = [1.0604174933173303e6, 7.464449633034707e6]

predict(PARS[end], p_fixed)[4, :]

loss(PARS[end], p_fixed)

res.u

res.objective

res.stats