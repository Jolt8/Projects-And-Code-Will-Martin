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
using ComponentArrays

struct CellData
    volume::Float64
    centroid::Vec{3, Float64}
    material_id::Int
end

struct Connection
    cell_idx_a::Int
    cell_idx_b::Int
    area::Float64
    normal::Vec{3, Float64}
    distance::Float64
end

abstract type AbstractBC end

struct HeatBC <: AbstractBC
    type::Symbol 
    initial::Float64
end

abstract type AbstractPhysics end

struct HeatPhysics <: AbstractPhysics
    k::Float64 
    rho::Float64
    cp::Float64
    source_term::Float64 #volumetric heating
end

struct FVMMesh
    mesh_cells::Vector{CellData}
    connections::Vector{Connection}
end


struct MultiPhysicsBCs
    temp1::Vector{HeatBC}
    temp2::Vector{HeatBC} 
end

struct BoundarySystem
    boundary_map::MultiPhysicsBCs
    free_idxs::Vector{Int}
    dirichlet_idxs::Vector{Int}
end

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


function get_k_effective(k_a, k_b, alpha_a, alpha_b, beta)
    #k_modified_a = k_void + (alpha_a^3) * (k_a - k_void) + 1e-6 # 1e-6 prevents divide by zero
    #k_modified_b = k_void + (alpha_b^3) * (k_b - k_void) + 1e-6
    #we will probably use SIMP like above later 

    η = 0.5 #\eta
    alpha_projected_a = (tanh(beta*η) + tanh(beta*(alpha_a - η))) / (tanh(beta·η) + tanh(beta*(1 - η)))
    alpha_projected_b = (tanh(beta*η) + tanh(beta*(alpha_b - η))) / (tanh(beta·η) + tanh(beta*(1 - η)))

    #for now, we use RAMP
    #usually q is roughly 3 to 5
    k_void = 0.001
    q = 4

    k_modified_a = k_void + ((alpha_projected_a * (1.0 + q) * (k_a - k_void)) / (1.0 + q * alpha_projected_a))
    k_modified_b = k_void + ((alpha_projected_b * (1.0 + q) * (k_b - k_void)) / (1.0 + q * alpha_projected_b))
    
    # Harmonic mean for interface conductivity
    k_effective = 2 * k_modified_a * k_modified_b / (k_modified_a + k_modified_b)

    return k_effective
end

function build_fvm_geometry(grid, cell_sets)
    n_cells = getncells(grid)

    ref_shape = getrefshape(grid.cells[1])
    poly_interp = Lagrange{ref_shape, 1}() #1 = linear elements 2 = quadratic/curved edges
    facet_qr = FacetQuadratureRule{ref_shape}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation
    facet_values = FacetValues(facet_qr, poly_interp)
    cell_qr = QuadratureRule{ref_shape}(2)
    cell_values = CellValues(cell_qr, poly_interp)

    cells_data = Vector{CellData}(undef, n_cells)
    
    for cell in CellIterator(grid) #You have to use a CellIterator Here
        cell_id = cellid(cell)
        Ferrite.reinit!(cell_values, cell) #not using Ferrite. before reinit breaks this! (DO NOT USE reinit! ON ITS OWN)
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values))
        coords = getcoordinates(cell)
        cent = sum(coords) / length(coords)

        mat_id = 1

        for i in eachindex(cell_sets)
            if cell_id in cell_sets[i].cell_set_idxs
                mat_id = cell_sets[i].mat_id
            end
        end

        cells_data[cell_id] = CellData(vol, cent, mat_id)
    end

    top = ExclusiveTopology(grid)
    connections = Vector{Connection}()
    for i in 1:n_cells
        for face_idx in 1:nfacets(grid.cells[i])
            neighbor_info = top.face_face_neighbor[i, face_idx]

            if !isempty(neighbor_info)
                neighbor_idx = collect(neighbor_info[1].idx)[1]
                
                if i < neighbor_idx
                    coords = getcoordinates(grid, i)
                    Ferrite.reinit!(facet_values, coords, face_idx) #not using Ferrite. before reinit breaks this! (DO NOT USE reinit! ON ITS OWN)

                    area = sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values))

                    n_ref = getnormal(facet_values, 1) 

                    centroid_a = cells_data[i].centroid
                    centroid_b = cells_data[neighbor_idx].centroid
                    dist = norm(centroid_b - centroid_a)
                    #TODO: fix for non-orthogonal meshes

                    push!(connections, Connection(i, neighbor_idx, area, n_ref, dist))
                end
            end
        end
    end
    
    return FVMMesh(cells_data, connections)
end

function FVM_iter_f!(du, u, p, geo::FVMMesh, bc_sys::BoundarySystem, heat_phys::Vector{HeatPhysics}, beta::Float64, ax) 
    #you could add another vector of fluid_phys if needed like fluid_phys::Vector{FluidPhysics}
    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)
    
    du .= 0.0

    for conn in geo.connections
        idx_a = conn.cell_idx_a
        idx_b = conn.cell_idx_b

        alpha_a = p[idx_a]
        alpha_b = p[idx_b]
        
        mat_a = geo.mesh_cells[idx_a].material_id
        mat_b = geo.mesh_cells[idx_b].material_id

        k_a = heat_phys[mat_a].k
        k_b = heat_phys[mat_b].k

        k_effective = get_k_effective(k_a, k_b, alpha_a, alpha_b, beta)

        #temp 1 
        u_a = u.temp1[idx_a]
        u_b = u.temp1[idx_b]

        F_1 = numerical_flux(k_effective, u_a, u_b, conn.area, conn.distance)
        
        du.temp1[idx_a] -= F_1
        du.temp1[idx_b] += F_1

        #temp 2
        u_a_2 = u.temp2[idx_a]
        u_b_2 = u.temp2[idx_b]

        F_2 = numerical_flux(k_effective, u_a_2, u_b_2, conn.area, conn.distance)

        du.temp2[idx_a] -= F_2
        du.temp2[idx_b] += F_2
    end

    # Source and Capacity Loop
    for i in bc_sys.free_idxs
        vol = geo.mesh_cells[i].volume
        mat = geo.mesh_cells[i].material_id
        
        rho = heat_phys[mat].rho
        cp  = heat_phys[mat].cp

        S = heat_phys[mat].source_term

        du.temp1[i] += S
        du.temp2[i] += S
        
        cap = capacity(rho, cp, vol)
        du.temp1[i] /= cap
        du.temp2[i] /= cap
    end

    for i in bc_sys.dirichlet_idxs
        du.temp1[i] = 0.0
        du.temp2[i] = 0.0
    end
end


grid_dimensions = (3, 3, 3) # Smaller grid for testing AD speed
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

cell_half_dist = (right[1] / grid_dimensions[1]) / 2 
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))

heated_copper_physics = HeatPhysics(401.0, 8960.0, 385.0, 0)#80000.0) #note that if we created a new struct for copper_physics performance would die
steel_physics = HeatPhysics(30.0, 8000.0, 460.0, 0)

heat_phys_vec = [heated_copper_physics, steel_physics]

copper_cell_set_idxs = Set(getcellset(grid, "copper"))
steel_cell_set_idxs = Set(getcellset(grid, "steel"))

struct CellSet
    mat_id::Int
    cell_set_idxs::Set
end

copper_set = CellSet(1, copper_cell_set_idxs)
steel_set = CellSet(2, steel_cell_set_idxs)

cell_sets = [copper_set, steel_set]

free_idxs = Int[]
dirichlet_idxs = Int[]

function my_bc_mapper(cell_id)
    if cell_id in copper_cell_set_idxs
        bc_type_a = HeatBC(:Neumann, 500.0)
        bc_type_b = HeatBC(:Neumann, 700.0)
        push!(free_idxs, cell_id)
        return [bc_type_a, bc_type_b] #use :Dirichlet to fix temperature to initial in HeatBC
    elseif cell_id in steel_cell_set_idxs
        bc_type_a = HeatBC(:Neumann, 300.0)
        bc_type_b = HeatBC(:Neumann, 400.0)
        push!(free_idxs, cell_id)
        return [bc_type_a, bc_type_b]
    end
end

temp1_bcs = HeatBC[]
temp2_bcs = HeatBC[]

n_cells = length(grid.cells)
n_vars = 2

for cell_id in 1:n_cells
    bcs = my_bc_mapper(cell_id) # Returns vector [BC1, BC2]
    push!(temp1_bcs, bcs[1])
    push!(temp2_bcs, bcs[2])
end

boundary_map = MultiPhysicsBCs(temp1_bcs, temp2_bcs)

bc_sys = BoundarySystem(boundary_map, free_idxs, dirichlet_idxs)

fvm_geo = build_fvm_geometry(grid, cell_sets)

p_guess = 0.5 .+ 0.01 .* randn(n_cells)

u_proto = ComponentArray(temp1 = zeros(n_cells), temp2 = zeros(n_cells))

for cell_id in 1:n_cells
    u_proto.temp1[cell_id] = bc_sys.boundary_map.temp1[cell_id].initial
    u_proto.temp2[cell_id] = bc_sys.boundary_map.temp2[cell_id].initial
end

u_axes = getaxes(u_proto)[1] 
#we use u_proto and u_axes in the function because KrylovJL_GMRES complains when a component arrays is passed in

u0 = Vector{Float64}(u_proto)

f_closure = (du, u, p) -> FVM_iter_f!(du, u, p, fvm_geo, bc_sys, heat_phys_vec, 1.0, u_axes)

detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
#not sure if pure TracerSparsityDetector is faster

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

# Implicit Solving Stuff Below (just add t to FVM_iter_f!)
#=
f_closure_implicit = (du, u, p, t) -> FVM_iter_f!(du, u, p, t, fvm_geo, bc_sys, heat_phys_vec, 1.0, u_axes)
#just re-add t to the FVM_iter_f! function above to make it compatible with implicit solving

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure_implicit(du, u, p_guess, 0.0), du0, u0, detector)

ode_func = ODEFunction(f_closure_implicit, jac_prototype = float.(jac_sparsity))

t0, tMax = 0.0, 1000.0 
desired_steps = 10
dt = tMax / desired_steps
tspan = (t0, tMax)
t = t0:dt:tMax;

implicit_prob = ODEProblem(ode_func, u0, tspan, p_guess)

function iluzero(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu0(convert(AbstractMatrix, W))
    else
        Pl = Plprev
    end
    Pl, nothing
end

@time sol = solve(implicit_prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true))
=#

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

cb = function (state, l)
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

reuse_prev_optimized_parameters = false

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
#Logging.disable_logging(Logging.Warn - 1)  # enable all warnings

for i in eachindex(beta_vals)
    f_closure = (du, u, p) -> FVM_iter_f!(du, u, p, fvm_geo, bc_sys, heat_phys_vec, 1.0, u_axes)
    
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
        return Array(sol)#, SciMLBase.successful_retcode(sol)
    end

    function loss(θ)
        #pred_u, is_successful = predict(θ) #if we ever need to return a high value for crashing the solver
        pred_u = predict(θ)

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
