using Ferrite
using DifferentialEquations
using LinearAlgebra
using Dates
using SparseArrays
using SciMLSensitivity
using Optimization, OptimizationPolyalgorithms, Zygote
using Enzyme
import SparseConnectivityTracer, ADTypes 

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

grid_dimensions = (6, 3, 3) # Smaller grid for testing AD speed
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
            boundary_map[i] = HeatBC(:Neumann, 1000.0)
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

function FVM_iter_f!(du, u, p, t, geo::FVMGeometry)
    du .= 0.0
    
    # Unpack parameters for clarity (Zygote handles index access well)
    # Material 1 (Copper)
    k1, rho1, cp1 = p[1], p[2], p[3]
    # Material 2 (Steel)
    k2, rho2, cp2 = p[4], p[5], p[6]

    # Internal Flux Loop
    # Note: iterating over struct arrays (geo.connections) is fine if they are constant
    for conn in geo.connections
        idx_a = conn.cell_idx_a
        idx_b = conn.cell_idx_b
        
        mat_a = geo.mesh_cells[idx_a].material_id
        mat_b = geo.mesh_cells[idx_b].material_id
        
        # Get k based on material ID
        k_a = (mat_a == 1) ? k1 : k2
        k_b = (mat_b == 1) ? k1 : k2
        
        # Harmonic mean for interface conductivity
        k_avg = 2 * k_a * k_b / (k_a + k_b)

        u_a = u[1, idx_a]
        u_b = u[1, idx_b]
        
        # Simplified flux call (removed normal for scalar heat eq unless anisotropic)
        F = numerical_flux(k_avg, u_a, u_b, conn.area, conn.distance)
        
        du[1, idx_a] -= F
        du[1, idx_b] += F
    end

    # Source and Capacity Loop
    for i in geo.free_idxs
        vol = geo.mesh_cells[i].volume
        mat = geo.mesh_cells[i].material_id
        
        # Get rho, cp based on material
        rho = (mat == 1) ? rho1 : rho2
        cp  = (mat == 1) ? cp1 : cp2
        
        # Assuming 0 source term for now, or add to p
        S = 0.0 

        du[1, i] += S
        
        cap = capacity(rho, cp, vol)
        du[1, i] /= cap
    end

    for i in geo.dirichlet_idxs
        du[1, i] = 0.0
    end
end

# Create the ODEProblem using a closure to capture `fvm_geo`
# The syntax is f(du, u, p, t)
f_closure = (du, u, p, t) -> FVM_iter_f!(du, u, p, t, fvm_geo)

n_cells = length(fvm_geo.mesh_cells)
u0 = zeros(1, n_cells)
for i in 1:n_cells
    u0[1, i] = fvm_geo.boundary_map[i].initial
end


t0, tMax = 0.0, 1000.0 
desired_steps = 100
dt = tMax / 100
tspan = (t0, tMax)
t = t0:dt:tMax;
#prob = ODEProblem(f_closure, u0, tspan, p_true)

println("Solving Forward Problem...")
# Using Tsit5 for explicit, ensure dt is small enough for stability
@time sol = solve(prob, Tsit5(), saveat=t)
stats = @timed sol = solve(prob, Tsit5(), saveat=t)
solve_time = stats.time
@time sol = solve(prob, Tsit5(), saveat=t)
arr_sol = Array(sol) # The "Truth" data to train against

#formatted arr_sol[1 (for some reason), cell_index, time_index]


# Start with wrong parameters

p_guess = log.([401.0, 30.0])
p_fixed = log.([8960.0, 385.0, 8000.0, 460.0])

function assemble_p(θ, constants)
    k1, k2 = θ[1], θ[2]
    rho1, cp1, rho2, cp2 = constants[1], constants[2], constants[3], constants[4]
    
    # Reconstruct the order expected by FVM_iter_f!
    # [k1, rho1, cp1, k2, rho2, cp2]
    return [k1, rho1, cp1, k2, rho2, cp2]
end

exp.(assemble_p(p_guess, p_fixed))

function predict(θ, constants)
    # solve needs to know to differentiate through θ (which enters as p)
    # sensealg=InterpolatingAdjoint() is standard for memory efficiency
    p_full = exp.(assemble_p(θ, constants))
    Array(solve(prob, Tsit5(), p=p_full, saveat=t, sensealg=InterpolatingAdjoint(autodiff=AutoEnzyme())))
end
#for some reason interpolating adjoint seems much faster than GaussAdjoint for some reason

#=
function loss(θ)
    pred = predict(θ)
    # Simple MSE
    return sum(abs2, pred .- (zeros(Float64, length(arr_sol[:, 1])) .+ 400.0))
end=#


function loss(θ)
    pred = predict(θ)

    #if length(pred[1, :]) < 101
        #return 1e9 # to prevent it from trying to crash the simulation and get its initial value
    #end

    vol_fraction = mean(θ)
    vol_target = 0.4
    vol_penalty = 10000.0 * max(0.0, vol_fraction - vol_target)^2
    #this vol_penalty thing never ever works
    
    #should probably check if its within the desired temperature for the entire simulation time 
    final_temp = pred[grid_idx_closest_to_desired_temp, end]
    physics_penalty = abs2(final_temp - desired_temp)
    
    K = (3529 * 0.5) #+ 3000.0 #it has half the effect on the loss than the physics themselves with the * 0.5
    #the 3127 is just the initial physics error
    #(physics_penalty + vol_penalty) * 100 
    #I think this will make it so that it always tries to get the optimimal solution first and then progress towards only 0 or 1 for alpha
    #don't do this, it just minimizes the volume and physics because K is depenent on those and if both of those = 0, alpha errors = 0
    #perhaps using K * the initial physics penalty would work

    #not entirely sure if the if θ[i]... is a good idea
    alpha_errors = sum(@. abs(16 * K * (θ^2 * (θ - 1)^2)))

    average_alpha_penalty = alpha_errors / length(θ)
    #not sure if averaging it by dividing by the length of θ is a good idea

    #total_penalty = physics_penalty + vol_penalty + alpha_errors
    total_penalty = physics_penalty + average_alpha_penalty #+ average_volume

    #println("physics_penalty, ", physics_penalty)
    #println("alpha_errors, ", average_alpha_penalty)

    return total_penalty
end

half_point = Int((grid_dimensions[1] * grid_dimensions[2] * grid_dimensions[3]) / 2)

test_pred = predict(p_guess, p_fixed)

test_pred[1, half_point, 50]

desired_temp = 900.0

sum(abs2, test_pred[1, half_point, 50] - desired_temp)

function loss(θ, constants)
    pred = predict(θ, constants)
    # Simple MSE
    #return (sum(abs2, pred[4, :] .- 400.0) / length(test_pred))
    return sum(abs2, pred[1, half_point, 50] - desired_temp) #desired temp at that point
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

println("Starting Optimization...")

# f_abstol: Stop if the Loss function changes less than this amount
# g_abstol: Stop if the Gradient is smaller than this (the slope is flat)

desired_improvement = 2.00 #x times better

corresponding_f_abstol = loss(p_guess, p_fixed) / desired_improvement 

#we do this just to ensure that it doesn't try to optimize further by forcing an impossible conductivity that makes the solver take forever 
#and take very, very, small time steps

res = Optimization.solve(optprob, PolyOpt(), callback=cb, 
                         f_abstol = corresponding_f_abstol, 
                         g_abstol = 1e-3,
)


res.stats

optimized_parameters_log = res.u

optimized_parameters = exp.(res.u)

predict(optimized_parameters, p_fixed)[1, half_point, :]

predict(optimized_parameters, p_fixed)[1, half_point, Int(desired_steps/2)] #it always optimizes for it to fail solving, ughh

loss(optimized_parameters, p_fixed)

res.u

res.objective
