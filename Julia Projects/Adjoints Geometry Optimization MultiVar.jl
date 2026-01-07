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

struct PhysicsBC

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

grid_dimensions = (3, 3, 3) # Smaller grid for testing AD speed
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

# --- 5. The ODE Function (Closure) ---
# We pass fvm_geo as a captured variable (closure). 
# It is NOT in 'p', so AD ignores it (which is good).

function FVM_iter_f!(du, u, p, geo::FVMGeometry, phys::PhysicsParams, beta::Float64, n_cells::Int)
    N = n_cells
    temp_2_idx = 1 * N
    #you could also do pressure_idx = 2 * N etc. 
    
    du .= 0.0

    #constants are located at 2nd index
    k1, rho1, cp1 = phys.k_cu, phys.rho_cu, phys.cp_cu
    k2, rho2, cp2 = phys.k_st, phys.rho_st, phys.cp_st

    # Internal Flux Loop
    # Note: iterating over struct arrays (geo.connections) is fine if they are constant
    for conn in geo.connections
        idx_a = conn.cell_idx_a
        idx_b = conn.cell_idx_b

        alpha_a = p[idx_a]
        alpha_b = p[idx_b]
        
        mat_a = geo.mesh_cells[idx_a].material_id
        mat_b = geo.mesh_cells[idx_b].material_id
        
        # Get k based on material ID
        k_a = (mat_a == 1) ? k1 : k2
        k_b = (mat_b == 1) ? k1 : k2

        #k_modified_a = k_void + (alpha_a^3) * (k_a - k_void) + 1e-6 # 1e-6 prevents divide by zero
        #k_modified_b = k_void + (alpha_b^3) * (k_b - k_void) + 1e-6
        #we will probably use SIMP like above later 

        η = 0.5 #\eta
        alpha_projected_a = (tanh(beta*η) + tanh(beta*(alpha_a - η))) / (tanh(beta·η) + tanh(beta*(1 - η)))
        alpha_projected_b = (tanh(beta*η) + tanh(beta*(alpha_b - η))) / (tanh(beta·η) + tanh(beta*(1 - η)))

        #1 → 1.5 → 2 → 3 → 4 → 6 → 8 → 12 → 16 → 24 → 32

        #for now, we use RAMP
        #usually q is roughly 3 to 5
        k_void = 0.001
        q = 4

        k_modified_a = k_void + ((alpha_projected_a * (1.0 + q) * (k_a - k_void)) / (1.0 + q * alpha_projected_a))
        k_modified_b = k_void + ((alpha_projected_b * (1.0 + q) * (k_b - k_void)) / (1.0 + q * alpha_projected_b))
        
        # Harmonic mean for interface conductivity
        k_avg = 2 * k_modified_a * k_modified_b / (k_modified_a + k_modified_b)

        u_a = u[idx_a]
        u_b = u[idx_b]
        u_a_2 = u[idx_a + temp_2_idx]
        u_b_2 = u[idx_b + temp_2_idx]
        #println("u_a_1, ", u_a, "idx: ", idx_a)
        #println("u_b_1, ", u_b, "idx: ", idx_b)
        #println("u_a_2, ", u_a_2, "idx: ", idx_a + 1 * n_cells)
        #println("u_b_2, ", u_b_2, "idx: ", idx_b + 1 * n_cells)
        
        # Simplified flux call (removed normal for scalar heat eq unless anisotropic)
        F_1 = numerical_flux(k_avg, u_a, u_b, conn.area, conn.distance)
        F_2 = numerical_flux(k_avg, u_a_2, u_b_2, conn.area, conn.distance)
        
        du[idx_a] -= F_1
        du[idx_b] += F_1
        du[idx_a + temp_2_idx] -= F_2
        du[idx_b + temp_2_idx] += F_2
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

        du[i] += S
        du[i + temp_2_idx] += S
        
        cap = capacity(rho, cp, vol)
        du[i] /= cap
        du[i + temp_2_idx] /= cap
    end

    for i in geo.dirichlet_idxs
        du[i] = 0.0
        du[i + temp_2_idx] = 0.0
    end
end

n_cells = length(grid.cells)
p_guess = 0.5 .+ 0.01 .* randn(n_cells)

phys_params = PhysicsParams(401.0, 8960.0, 385.0, 30.0, 8000.0, 460.0)

# Create the ODEProblem using a closure to capture `fvm_geo`
# The syntax is f(du, u, p, t)
length(fvm_geo.mesh_cells)
f_closure = (du, u, p) -> FVM_iter_f!(du, u, p, fvm_geo, phys_params, 1.0, length(fvm_geo.mesh_cells))

n_cells = length(fvm_geo.mesh_cells)
n_vars = 2
#u0 = zeros(n_cells * n_vars)

u0 = zeros(n_cells * n_vars)

for i in 1:n_cells
    for j in 1:n_vars
        println(i + (n_cells * (j - 1)))
        u0[i + (n_cells * (j - 1))] = fvm_geo.boundary_map[i].initial[j]
    end
end
u0

detector = SparseConnectivityTracer.TracerSparsityDetector()

du0 = u0 .* 0.0

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

n_cells = length(grid.cells)
p_guess = 0.5 .+ 0.01 .* randn(n_cells)

println("Solving Forward Problem...")
#@time sol = solve(prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), saveat=t)
@time sol = solve(prob, NonlinearSolve.NewtonRaphson(concrete_jac = true), p=p_guess)

stop = 1

#=
f_closure_implicit = (du, u, p, t) -> FVM_iter_f!(du, u, p, t, fvm_geo, phys_params, 1.0)
#adding t to make it work with the implicit solve to test stability

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0, tMax = 0.0, 1000.0 
desired_steps = 10
dt = tMax / desired_steps
tspan = (t0, tMax)
t = t0:dt:tMax;

implicit_prob = ODEProblem(ode_func, u0, tspan, p_guess)

#@time sol = solve(prob, Tsit5())

@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), saveat=t)
=#



#=
function test_predict(θ)
    #sol = solve(prob, Tsit5(), p=θ, saveat=t, sensealg=InterpolatingAdjoint(autodiff=AutoEnzyme()))
    sol = solve(prob, NonlinearSolve.NewtonRaphson(), p=θ, sensealg=SteadyStateAdjoint(autodiff=AutoEnzyme()), verbose = false)
    return Array(sol), SciMLBase.successful_retcode(sol)
end

#for some reason interpolating adjoint seems much faster than GaussAdjoint for some reason 
#even though interpolatingAdjoint is reported to be more computationally expensive

#use SteadyStateAdjoint here

#SteadyStateProblem with DynamicSS: Good if your system might have stability issues or if you want to leverage your existing ODE infrastructure. 
#Also good if you might later want transient optimization.

#SteadyStateProblem with SSRootfind or NonlinearProblem: Better for pure steady-state. Faster when it works. 
#Requires a good initial guess and a well-conditioned Jacobian.

test_pred = test_predict(p_guess)[1]

desired_temp = 900.0

grid_idx_closest_to_desired_temp = argmin(abs.(test_pred .- desired_temp))

#test_pred[grid_idx_closest_to_desired_temp, end]

function physics_loss(pred_array)
    final_temp = pred_array[grid_idx_closest_to_desired_temp]
    physics_penalty = abs2(final_temp - desired_temp)
    return physics_penalty
end

initial_physics_loss = physics_loss(test_pred)

function alpha_loss(θ)
    K = (3529 * 0.5)
    alpha_errors = sum(@. abs(16 * K * (θ^2 * (θ - 1)^2)))

    average_alpha_penalty = alpha_errors / length(θ)
    return average_alpha_penalty
end

initial_alpha_loss = alpha_loss(p_guess)

#test_pred = predict(p_guess)[1]
#initial_physics_loss = physics_loss(test_pred)

function test_loss(θ)
    pred_u, is_successful = test_predict(θ)

    breaking_solver_penalty = 0.0

    if is_successful != true
        println("hoo dog, we messed up")
        breaking_solver_penalty = 1.0e9
    end
    
    scaled_physics_loss = physics_loss(pred_u) / initial_physics_loss 

    scaled_alpha_loss = alpha_loss(θ) / initial_alpha_loss

    #println("scaled physics loss, ", scaled_physics_loss)

    #println("scaled alpha loss, ", scaled_alpha_loss)

    return (1.0 * scaled_physics_loss) + (1.0 * scaled_alpha_loss) #+ breaking_solver_penalty
end


#l = loss(p_guess)

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
    verbose = false
end

#adtype = Optimization.AutoZygote()
#optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
#optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype, cons_jac_prototype=jac_sparsity)
#adding jac_sparsity kills performance
#update: it REALLY kills performance, I'm probably doing something wrong somewhere
#update 2: apparently passing the jac_prototype here is actually uncessary and the only place in which you need to pass it is in the original implicit ODE solve

#optprob = Optimization.OptimizationProblem(optf, p_guess, lb=lower_bounds, ub=upper_bounds, verbose = false)
#println("Starting Optimization...")

# f_abstol: Stop if the Loss function changes less than this amount
# g_abstol: Stop if the Gradient is smaller than this (the slope is flat)

SciMLSensitivity.STACKTRACE_WITH_VJPWARN[] = false #turn to true to debug enzyme JVP to work

lower_bounds = zeros(n_cells)
upper_bounds = ones(n_cells)

beta_steps = 12

beta_vals = []

for i in 1:beta_steps
    if i == 1
        push!(beta_vals, 1.0)
    else
        push!(beta_vals, 2 * beta_vals[i-1])
    end
end

beta_vals

#beta_vals = [1.0, 2.0]

seeing_if_betas_are_used = []

reuse_prev_optimized_parameters = true

#kinda pointless, but whatever
if reuse_prev_optimized_parameters 
    initial_parameters = new_guess_params
else 
    initial_parameters = p_guess
end

new_guess_params = 0


#ode_prob = ODEProblem(f_closure, u0, tspan, p_guess)

#prob = NonlinearProblem(ode_prob)

Logging.disable_logging(Logging.Warn)  # Disable all warnings
#Logging.disable_logging(Logging.Warn - 1)  # Disable all warnings

for i in eachindex(beta_vals)
    f_closure = (du, u, p) -> FVM_iter_f!(du, u, p, fvm_geo, phys_params, beta_vals[i], length(fvm_geo.mesh_cells))
    
    if i == 1 
        detector = SparseConnectivityTracer.TracerSparsityDetector()
        jac_sparsity = ADTypes.jacobian_sparsity(
            (du, u) -> f_closure(du, u, p_guess), du0, u0, detector)
        seeing_if_betas_are_used = []
    end

    println(beta_vals[i])
    push!(seeing_if_betas_are_used, beta_vals[i])

    #currently for a 3x3x3 system adding the jac_prototype makes it 3x slower but that will probably change for larger systems
    
    nl_func = NonlinearFunction(f_closure, jac_prototype = float.(jac_sparsity))

    prob = NonlinearProblem(nl_func, u0, p_guess)

    function predict(θ)
        sol = solve(prob, NonlinearSolve.NewtonRaphson(), p=θ, sensealg=SteadyStateAdjoint(diff_type=AutoEnzyme()), verbose = false, alias=SciMLBase.NonlinearAliasSpecifier(alias_u0=true))
        return Array(sol), SciMLBase.successful_retcode(sol)
    end

    function loss(θ)
        pred_u, is_successful = predict(θ)

        scaled_physics_loss = physics_loss(pred_u) / initial_physics_loss 

        scaled_alpha_loss = alpha_loss(θ) / initial_alpha_loss
        println("scaled_physics_loss, ", scaled_physics_loss)
        println("scaled_alpha_loss, ", scaled_alpha_loss)

        #println("scaled physics loss, ", scaled_physics_loss)

        #println("scaled alpha loss, ", scaled_alpha_loss)

        return (1.0 * scaled_physics_loss) + (1.0 * scaled_alpha_loss) #+ breaking_solver_penalty
        #for some reason, the only way to make it optimize towards 0 and 1 is to add the scaled alpha loss.
        #the tanh()... for alpha doesn't seem to do this reliably
    end

    if i == 1
        new_guess_params = initial_parameters
    end

    adtype = Optimization.AutoZygote()
    optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
    #optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype, cons_jac_prototype=jac_sparsity)
    #adding jac_sparsity kills performance
    #update 2: apparently passing the jac_prototype here is actually uncessary and the only place in which you need to pass it is in the original implicit ODE solve

    lower_bounds = zeros(n_cells)
    upper_bounds = ones(n_cells)

    optprob = Optimization.OptimizationProblem(optf, new_guess_params, lb=lower_bounds, ub=upper_bounds, verbose=false)

    @time res = Optimization.solve(optprob, OptimizationOptimJL.Fminbox(), 
                            callback=cb,
                            f_abstol = 1e-17, 
                            g_abstol = 1e-17,
                            verbose=false
    )
    
    println(res.u)

    new_guess_params = res.u
end

seeing_if_betas_are_used

optimized_parameters = new_guess_params

#VERY IMPORTANT!! Previously, u was formatted [1, cell_idx] because I wanted to be able to handle multiple u values per cell but 
#trying to get and sort of JVP algorithim to work with the solver was impossible because it complained that [1, 27] did not match the JVP size of [27, 27]

mean(optimized_parameters)

test_predict(optimized_parameters)[1]

test_predict(p_guess)[1][grid_idx_closest_to_desired_temp, end]

test_predict(optimized_parameters)[1][grid_idx_closest_to_desired_temp, end]

test_loss(p_guess)

test_loss(optimized_parameters)
=#
