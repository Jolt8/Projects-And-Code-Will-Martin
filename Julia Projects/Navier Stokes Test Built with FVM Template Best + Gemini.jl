using Ferrite
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using SciMLSensitivity
using Optimization, OptimizationPolyalgorithms
using Enzyme
using RecursiveArrayTools
using OptimizationOptimJL
using ILUZero
import AlgebraicMultigrid
import SparseConnectivityTracer, ADTypes 
using NonlinearSolve
import Logging
using ComponentArrays

# --- GEOMETRY STRUCTS (Unchanged) ---
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

struct FVMMesh
    mesh_cells::Vector{CellData}
    connections::Vector{Connection}
end

# --- PHYSICS & BC STRUCTS (Updated for Fluids) ---

abstract type AbstractBC end

struct FluidBC <: AbstractBC
    type::Symbol   # :Dirichlet or :Neumann (usually pressure outlet)
    value::Float64
end

struct FluidPhysics
    mu::Float64   # Dynamic Viscosity
    rho::Float64  # Density
    alpha_max::Float64 # Maximum porous resistance (Brinkman penalty for solids)
end

struct FluidSystemBCs
    vx::Vector{Union{FluidBC, Nothing}}
    vy::Vector{Union{FluidBC, Nothing}}
    vz::Vector{Union{FluidBC, Nothing}}
    p::Vector{Union{FluidBC, Nothing}}
end

struct BoundarySystem
    boundary_map::FluidSystemBCs
    # Indices where we solve (Free) vs Fixed
    dirichlet_idxs_vx::Vector{Int}
    dirichlet_idxs_vy::Vector{Int}
    dirichlet_idxs_vz::Vector{Int}
    dirichlet_idxs_p::Vector{Int}
end

# --- HELPER FUNCTIONS ---

# Upwind scheme for Convection
# If flow goes A -> B, carry property of A. Else B.
function upwind_val(u_a, u_b, mass_flux)
    return mass_flux > 0 ? u_a : u_b
end

# Effective Viscosity / Penalization for Topology Optimization
# alpha (design variable) = 0 (Liquid), 1 (Solid)
function get_material_props(phys::FluidPhysics, design_param, beta)
    # Interpolation (RAMP or SIMP)
    # 0 = Fluid, 1 = Solid
    # In fluid, resistance is 0. In solid, resistance is high.
    
    # Simple penalization: alpha_brinkman = alpha_max * q * p / (1 + q * p) 
    # Here we assume design_param is 0 for fluid, 1 for solid
    
    q = 10.0 # Tuning parameter for RAMP
    brinkman_alpha = phys.alpha_max * (q * design_param) / (1.0 + q * design_param)
    
    return phys.mu, phys.rho, brinkman_alpha
end

# --- GEOMETRY BUILDER (Unchanged) ---
function build_fvm_geometry(grid)
    n_cells = getncells(grid)
    ref_shape = getrefshape(grid.cells[1])
    poly_interp = Lagrange{ref_shape, 1}()
    facet_qr = FacetQuadratureRule{ref_shape}(2)
    facet_values = FacetValues(facet_qr, poly_interp)
    cell_qr = QuadratureRule{ref_shape}(2)
    cell_values = CellValues(cell_qr, poly_interp)

    cells_data = Vector{CellData}(undef, n_cells)
    
    for cell in CellIterator(grid)
        cell_id = cellid(cell)
        Ferrite.reinit!(cell_values, cell)
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values))
        coords = getcoordinates(cell)
        cent = sum(coords) / length(coords)
        cells_data[cell_id] = CellData(vol, cent, 1) # Single material ID for now
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
                    Ferrite.reinit!(facet_values, coords, face_idx)
                    area = sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values))
                    n_ref = getnormal(facet_values, 1)
                    dist = norm(cells_data[neighbor_idx].centroid - cells_data[i].centroid)
                    push!(connections, Connection(i, neighbor_idx, area, n_ref, dist))
                end
            end
        end
    end
    return FVMMesh(cells_data, connections)
end

# --- NAVIER STOKES RESIDUAL FUNCTION ---

function NS_iter_f!(du, u, p_topo, geo::FVMMesh, bc_sys::BoundarySystem, phys::FluidPhysics, ax) 
    # Unpack
    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)
    
    # Reset residuals
    du .= 0.0

    # 1. Flux Loop (Interior Faces)
    for conn in geo.connections
        idx_a, idx_b = conn.cell_idx_a, conn.cell_idx_b
        dist, area, n = conn.distance, conn.area, conn.normal

        # Interpolate variables to face (Linear)
        vx_f = 0.5 * (u.vx[idx_a] + u.vx[idx_b])
        vy_f = 0.5 * (u.vy[idx_a] + u.vy[idx_b])
        vz_f = 0.5 * (u.vz[idx_a] + u.vz[idx_b])
        
        # --- CONTINUITY & PRESSURE STABILIZATION ---
        # Mass Flux (mdot) = rho * area * (V . n)
        v_dot_n = vx_f*n[1] + vy_f*n[2] + vz_f*n[3]
        mass_flux = phys.rho * area * v_dot_n

        # Rhie-Chow / Stabilization:
        # Standard colocated FVM requires adding a term proportional to (P_a - P_b) 
        # to the mass flux to prevent checkerboarding.
        # Stabil_factor ~ rho * Area / coeff_momentum. Simplified here:
        epsilon = 1e-4 * (area / dist) # Tuning parameter
        p_diff = u.p[idx_b] - u.p[idx_a]
        stabilization_flux = -epsilon * p_diff 

        total_mass_flux = mass_flux + stabilization_flux

        # Add to Continuity Residual (Conservation of Mass)
        # Note: If mass leaves A, residual decreases? 
        # Standard: Div(u) = 0. Outward flux is positive.
        du.p[idx_a] += total_mass_flux
        du.p[idx_b] -= total_mass_flux

        # --- MOMENTUM EQUATION ---
        
        # 1. Diffusion (Viscous Stress)
        # mu * Area * dU/dx
        diff_factor = phys.mu * area / dist
        
        diff_vx = diff_factor * (u.vx[idx_b] - u.vx[idx_a])
        diff_vy = diff_factor * (u.vy[idx_b] - u.vy[idx_a])
        diff_vz = diff_factor * (u.vz[idx_b] - u.vz[idx_a])

        du.vx[idx_a] += diff_vx
        du.vx[idx_b] -= diff_vx
        du.vy[idx_a] += diff_vy
        du.vy[idx_b] -= diff_vy
        du.vz[idx_a] += diff_vz
        du.vz[idx_b] -= diff_vz

        # 2. Convection (rho * u * u)
        # Flux = MassFlux * UpwindVelocity
        conv_vx = total_mass_flux * upwind_val(u.vx[idx_a], u.vx[idx_b], total_mass_flux)
        conv_vy = total_mass_flux * upwind_val(u.vy[idx_a], u.vy[idx_b], total_mass_flux)
        conv_vz = total_mass_flux * upwind_val(u.vz[idx_a], u.vz[idx_b], total_mass_flux)

        du.vx[idx_a] -= conv_vx # Leaving A
        du.vx[idx_b] += conv_vx # Entering B
        du.vy[idx_a] -= conv_vy
        du.vy[idx_b] += conv_vy
        du.vz[idx_a] -= conv_vz
        du.vz[idx_b] += conv_vz

        # 3. Pressure Gradient Force
        # Force = Pressure * Area * Normal
        # Pressure is scalar, acts normal to face
        p_avg = 0.5 * (u.p[idx_a] + u.p[idx_b])
        pressure_force_x = p_avg * area * n[1]
        pressure_force_y = p_avg * area * n[2]
        pressure_force_z = p_avg * area * n[3]

        du.vx[idx_a] -= pressure_force_x
        du.vx[idx_b] += pressure_force_x
        du.vy[idx_a] -= pressure_force_y
        du.vy[idx_b] += pressure_force_y
        du.vz[idx_a] -= pressure_force_z
        du.vz[idx_b] += pressure_force_z
    end

    # 2. Volumetric Loop (Source Terms & Penalization)
    for i in 1:length(geo.mesh_cells)
        vol = geo.mesh_cells[i].volume
        
        # Brinkman Penalization for Topology Optimization
        # S = - alpha * u
        # If p_topo[i] is 1 (solid), alpha is high -> u forced to 0
        mu, rho, alpha_brinkman = get_material_props(phys, p_topo[i], 1.0)
        
        damping = alpha_brinkman * vol

        du.vx[i] -= damping * u.vx[i]
        du.vy[i] -= damping * u.vy[i]
        du.vz[i] -= damping * u.vz[i]
    end

    # 3. Apply Boundary Conditions (Dirichlet)
    # The residuals for Dirichlet nodes are replaced by (u_i - BC_val)
    
    for idx in bc_sys.dirichlet_idxs_vx
        du.vx[idx] = u.vx[idx] - bc_sys.boundary_map.vx[idx].value
    end
    for idx in bc_sys.dirichlet_idxs_vy
        du.vy[idx] = u.vy[idx] - bc_sys.boundary_map.vy[idx].value
    end
    for idx in bc_sys.dirichlet_idxs_vz
        du.vz[idx] = u.vz[idx] - bc_sys.boundary_map.vz[idx].value
    end
    for idx in bc_sys.dirichlet_idxs_p
        du.p[idx]  = u.p[idx]  - bc_sys.boundary_map.p[idx].value
    end
end

# --- SETUP SIMULATION ---

# 1. Grid (Lid Driven Cavity Dimensions)
grid_dim = (10, 10, 10) # Keep small for direct solver speed
left = Vec{3}((0.0, 0.0, 0.0))
right = Vec{3}((1.0, 1.0, 0.2)) # Thin slice 3D
grid = generate_grid(Hexahedron, grid_dim, left, right)

# 2. Physics (Water-like)
# alpha_max = 1e6 (High resistance for solids)
fluid_phys = FluidPhysics(0.01, 1.0, 1e6) 

# 3. Boundary Conditions Logic
# Lid Driven Cavity: 
# Top (y=1.0) moves at Vx=1.0
# Others are No-Slip (V=0)
# Fix Pressure at one corner to 0 (make system determined)

n_cells = getncells(grid)
vx_bcs = Vector{Union{FluidBC, Nothing}}(nothing, n_cells)
vy_bcs = Vector{Union{FluidBC, Nothing}}(nothing, n_cells)
vz_bcs = Vector{Union{FluidBC, Nothing}}(nothing, n_cells)
p_bcs  = Vector{Union{FluidBC, Nothing}}(nothing, n_cells)

dir_vx, dir_vy, dir_vz, dir_p = Int[], Int[], Int[], Int[]

dx = (right[1] - left[1]) / grid_dim[1]
dy = (right[2] - left[2]) / grid_dim[2]
dz = (right[3] - left[3]) / grid_dim[3]

# Identify boundary cells via Ferrite (simplification: simple geometry check)
for cell in CellIterator(grid)
    idx = cellid(cell)
    center = sum(getcoordinates(cell))/8
    
    # 1. LID DRIVEN TOP (Moving Wall)
    # Check if we are within half a cell width of the top
    if center[2] > (right[2] - dy * 0.6)
        vx_bcs[idx] = FluidBC(:Dirichlet, 1.0) # VELOCITY IS 1.0 HERE
        vy_bcs[idx] = FluidBC(:Dirichlet, 0.0)
        vz_bcs[idx] = FluidBC(:Dirichlet, 0.0)
        push!(dir_vx, idx); push!(dir_vy, idx); push!(dir_vz, idx)
        
    # 2. OTHER WALLS (No Slip)
    elseif center[2] < (left[2] + dy * 0.6) || 
           center[1] < (left[1] + dx * 0.6) || 
           center[1] > (right[1] - dx * 0.6) ||
           center[3] < (left[3] + dz * 0.6) || 
           center[3] > (right[3] - dz * 0.6)
           
        vx_bcs[idx] = FluidBC(:Dirichlet, 0.0)
        vy_bcs[idx] = FluidBC(:Dirichlet, 0.0)
        vz_bcs[idx] = FluidBC(:Dirichlet, 0.0)
        push!(dir_vx, idx); push!(dir_vy, idx); push!(dir_vz, idx)
    end

    # 3. PRESSURE PIN (Pin one corner to 0 to prevent drift)
    if idx == 1
        p_bcs[idx] = FluidBC(:Dirichlet, 0.0)
        push!(dir_p, idx)
    end
end

boundary_sys = BoundarySystem(
    FluidSystemBCs(vx_bcs, vy_bcs, vz_bcs, p_bcs),
    dir_vx, dir_vy, dir_vz, dir_p
)

geo = build_fvm_geometry(grid)

# --- SOLVE ---

# Design Variable (Topology)
# 0 = Fluid, 1 = Solid. Initialize as all fluid.
p_topo = zeros(n_cells) 

# Initial Guess (Field Variables)
u_proto = ComponentArray(
    vx = zeros(n_cells), 
    vy = zeros(n_cells), 
    vz = zeros(n_cells), 
    p  = zeros(n_cells)
)

u0 = Vector{Float64}(u_proto)
count(!iszero, u0)
u_ax = getaxes(u_proto)[1]

# Closure
f_ns = (du, u, p) -> NS_iter_f!(du, u, p_topo, geo, boundary_sys, fluid_phys, u_ax)

# Sparsity Detection (Crucial for Coupled N-S)
println("Detecting Sparsity Pattern...")
detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
du0 = u0 .* 0.0
jac_sparsity = ADTypes.jacobian_sparsity((du, u) -> f_ns(du, u, p_topo), du0, u0, detector)

println("Initializing Solver...")
nl_func = NonlinearFunction(f_ns, jac_prototype = float.(jac_sparsity))
prob = NonlinearProblem(nl_func, u0)

function algebraicmultigrid(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(convert(AbstractMatrix, W)))
    else
        Pl = Plprev
    end
    Pl, nothing
end

println("Solving Navier-Stokes (Newton-Raphson)...")
# NewtonRaphson is robust for this monolithic formulation
@time sol = solve(prob, NonlinearSolve.NewtonRaphson(linsolve = LinearSolve.KrylovJL_GMRES(), concrete_jac = true), abstol=1e-5, reltol=1e-5)

# Extract results
res_u = ComponentVector(sol.u, u_ax)
println("Max Velocity X: ", maximum(res_u.vx))
println("Max Pressure: ", maximum(res_u.p))
println("Min Pressure: ", minimum(res_u.p))
record_sol = false

length(sol.u[1, :])

using WriteVTK

if record_sol == true
    date_and_time = Dates.format(now(), "I.MM.SS p yyyy-mm-dd")
    #date_and_time = Dates.format(now(), "I.MM.SS p")

    root_dir = "C://Users//wille//Desktop//Julia_cfd_output_files"

    project_name = replace(basename(@__FILE__),r".jl" => "")

    sim_folder_name = project_name * " " * date_and_time

    output_dir = joinpath(root_dir, sim_folder_name)

    mkpath(output_dir)

    pvd_filename = joinpath(output_dir, "solution_collection")

    pvd = paraview_collection(pvd_filename)

    step_filename = joinpath(output_dir, "timestep")

    #this step may actually become a significant bottleneck
    #update: This is a very big bottleneck
    for (step, t) in enumerate(sol.t)
        temp_field = sol.u[step][1, :]
        VTKGridFile(step_filename * " $step" * " at $(t)s.vtu", grid) do vtk
            write_cell_data(vtk, temp_field, "T")
            pvd[t] = vtk
        end
    end
    vtk_save(pvd)
end