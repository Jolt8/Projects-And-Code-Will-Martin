using Ferrite
using DifferentialEquations
using LinearAlgebra 
using WriteVTK
using Dates
using SparseArrays
using RecursiveArrayTools
using StaticArrays
import SparseConnectivityTracer, ADTypes 
using ProfileView

# ==============================================================================
# 1. CORE FVM ARCHITECTURE
# ==============================================================================

abstract type AbstractPhysics end

struct CellData
    volume::Float64
    centroid::Vec{3, Float64}
end

struct Connection
    cell_idx_a::Int
    cell_idx_b::Int
    area::Float64
    normal::Vec{3, Float64}
    distance::Float64
end

abstract type AbstractBC end

struct FVMProblem{B<:AbstractBC}
    mesh_cells::Vector{CellData}
    connections::Vector{Connection}
    boundary_map::Vector{B}
end

# ==============================================================================
# 2. NAVIER-STOKES (EULER) PHYSICS DEFINITION
# ==============================================================================

struct CompressibleEulerPhysics <: AbstractPhysics
    gamma::Float64  # Heat capacity ratio (1.4 for air)
    R_gas::Float64  # Specific gas constant
end

struct FluidBC <: AbstractBC
    type::Symbol # :Dirichlet (Inlet), :Neumann (Outlet), :Wall
    # State is [rho, rho*u, rho*v, rho*w, E]
    state::Vector{Float64} 
    physics::CompressibleEulerPhysics
end

# --- HELPER FUNCTIONS FOR FLUIDS ---

# Calculate Pressure from Conserved Variables
# p = (gamma - 1) * (E - 0.5 * rho * v^2)
function get_pressure(phys::CompressibleEulerPhysics, U)
    rho = U[1]
    rhou = U[2]
    rhov = U[3]
    rhow = U[4]
    E = U[5]
    
    # velocity squared * rho
    rho_v2 = (rhou^2 + rhov^2 + rhow^2) / rho
    
    return (phys.gamma - 1.0) * (E - 0.5 * rho_v2)
end

# --- THE RUSANOV FLUX FUNCTION ---
# This is the industry standard "Robust" flux for compressible flow.
# It stabilizes the solution by adding artificial dissipation at shocks.

function numerical_flux(phys::CompressibleEulerPhysics, u_L, u_R, area, normal, dist)
    T = eltype(u_L)
    # Unpack Left State
    rho_L = u_L[1]

    v_Lx = u_L[2]/rho_L
    v_Ly = u_L[3]/rho_L
    v_Lz = u_L[4]/rho_L

    p_L   = get_pressure(phys, u_L)
    E_L   = u_L[5]
    H_L   = (E_L + p_L) / rho_L # Enthalpy

    # Unpack Right State
    rho_R = u_R[1]

    v_Rx = u_R[2]/rho_R
    v_Ry = u_R[3]/rho_R
    v_Rz = u_R[4]/rho_R

    p_R   = get_pressure(phys, u_R)
    E_R   = u_R[5]
    H_R   = (E_R + p_R) / rho_R

    # Calculate Normal Velocities
    vn_L = v_Lx * normal[1] + v_Ly * normal[2] + v_Lz * normal[3]
    vn_R = v_Rx * normal[1] + v_Ry * normal[2] + v_Rz * normal[3]

    # 1. Compute Physical Fluxes (F_L and F_R)
    # Flux vector for Euler: [rho*vn, rho*u*vn + p*nx, ..., rho*H*vn]
    
    # Left Flux
    F_L = SVector{5, T}(
        rho_L * vn_L, 
        rho_L * v_Lx * vn_L + p_L * normal[1], 
        rho_L * v_Ly * vn_L + p_L * normal[2], 
        rho_L * v_Lz * vn_L + p_L * normal[3], 
        rho_L * H_L * vn_L
    ) 
    #F_L_1 = rho_L * vn_L
    #F_L_2 = rho_L * v_Lx * vn_L + p_L * normal[1]
    #F_L_3 = rho_L * v_Ly * vn_L + p_L * normal[2] 
    #F_L_4 = rho_L * v_Lz * vn_L + p_L * normal[3] 
    #F_L_5 = rho_L * H_L * vn_L
    #this may or may not be better
    
    # Right Flux
    F_R = SVector{5, T}(
        rho_R * vn_R, 
        rho_R * v_Rx * vn_R + p_R * normal[1], 
        rho_R * v_Ry * vn_R + p_R * normal[2], 
        rho_R * v_Rz * vn_R + p_R * normal[3], 
        rho_R * H_R * vn_R
    )
    #F_R_1 = rho_R * vn_R
    #F_R_2 = rho_R * v_Rx * vn_R + p_R * normal[1]
    #F_R_3 = rho_R * v_Ry * vn_R + p_R * normal[2]
    #F_R_4 = rho_R * v_Rz * vn_R + p_R * normal[3]
    #F_R_5 = rho_R * H_R * vn_R
    
    # 2. Rusanov Stabilization (The "Riemann Solver" part)
    # We need the maximum wave speed at the interface.
    # c = sqrt(gamma * p / rho)
    c_L = sqrt(abs(phys.gamma * p_L / rho_L))
    c_R = sqrt(abs(phys.gamma * p_R / rho_R))
    
    # Max wave speed
    lambda_max = max(abs(vn_L) + c_L, abs(vn_R) + c_R)

    # 3. Final Flux
    # F_face = 0.5 * (F_L + F_R) - 0.5 * lambda * (U_R - U_L)
    flux = SVector{5, T}(
        0.5 * (F_L[1] + F_R[1]) - 0.5 * lambda_max * (u_R[1] - u_L[1]), 
        0.5 * (F_L[2] + F_R[2]) - 0.5 * lambda_max * (u_R[2] - u_L[2]),
        0.5 * (F_L[3] + F_R[3]) - 0.5 * lambda_max * (u_R[3] - u_L[3]),
        0.5 * (F_L[4] + F_R[4]) - 0.5 * lambda_max * (u_R[4] - u_L[4]),
        0.5 * (F_L[5] + F_R[5]) - 0.5 * lambda_max * (u_R[5] - u_L[5])
    )
    #flux_1 = 0.5 * (F_L_1 + F_R_1) - 0.5 * lambda_max * (u_R[1] - u_L[1])
    #flux_2 = 0.5 * (F_L_2 + F_R_2) - 0.5 * lambda_max * (u_R[2] - u_L[2])
    #flux_3 = 0.5 * (F_L_3 + F_R_3) - 0.5 * lambda_max * (u_R[3] - u_L[3])
    #flux_4 = 0.5 * (F_L_4 + F_R_4) - 0.5 * lambda_max * (u_R[4] - u_L[4])
    #flux_5 = 0.5 * (F_L_5 + F_R_5) - 0.5 * lambda_max * (u_R[5] - u_L[5])
    #I tried doing this flux_1, flux_2, flux_3 etc. one but it didn't seem to improve the runtime at all

    return flux .* area #SVector{5, T}(flux_1 * area, flux_2 * area, flux_3 * area, flux_4 * area, flux_5 * area)
end

# No volumetric source for standard shock tube
function source(phys::CompressibleEulerPhysics, u, vol)
    return zeros(5) 
end

# For Euler, dU/dt = -Flux / Vol. So capacity is just volume.
function capacity(phys::CompressibleEulerPhysics, vol)
    return vol
end

function build_fvm_problem(grid, bc_map_func)
    n_cells = getncells(grid)
    
    # Ferrite Geometry Setup
    ref_shape = getrefshape(grid.cells[1])
    poly_interp = Lagrange{ref_shape, 1}() #1 = linear elements 2 = quadratic/curved edges
    facet_qr = FacetQuadratureRule{ref_shape}(2)#2 represents the number of integration points. Basically higher number = higher accuracy but more computation
    facet_values = FacetValues(facet_qr, poly_interp)
    cell_qr = QuadratureRule{ref_shape}(2)
    cell_values = CellValues(cell_qr, poly_interp)

    cells_data = Vector{CellData}(undef, n_cells)
    
    for cell in CellIterator(grid)
        id = cellid(cell)
        Ferrite.reinit!(cell_values, cell)
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values))
        coords = getcoordinates(cell)
        cent = sum(coords) / length(coords)
        cells_data[id] = CellData(vol, cent)
    end

    top = ExclusiveTopology(grid)
    connections = Vector{Connection}()
    boundary_map = Vector{FluidBC}(undef, n_cells)
    
    for i in 1:n_cells
        boundary_map[i] = bc_map_func(grid, i)

        for face_idx in 1:nfacets(grid.cells[i])
            neighbor_info = top.face_face_neighbor[i, face_idx]
            
            if !isempty(neighbor_info)
                neighbor_idx = collect(neighbor_info[1].idx)[1]
                
                if i < neighbor_idx
                    coords = getcoordinates(grid, i)
                    Ferrite.reinit!(facet_values, coords, face_idx)
                    area = sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values))
                    n_ref = getnormal(facet_values, 1) 
                    
                    cent_a = cells_data[i].centroid
                    cent_b = cells_data[neighbor_idx].centroid
                    dist = norm(cent_b - cent_a)

                    if (n_ref ⋅ (cent_b - cent_a)) < 0
                        n_ref = -n_ref
                    end

                    push!(connections, Connection(i, neighbor_idx, area, n_ref, dist))
                end
            end
        end
    end
    
    return FVMProblem(cells_data, connections, boundary_map)
end

function FVM_iter_f!(du, u, p::FVMProblem, t)
    du .= 0.0

    # Internal Faces
    for conn in p.connections
        idx_a = conn.cell_idx_a
        idx_b = conn.cell_idx_b
        
        # Views into Matrix columns
        u_a = @view u[:, idx_a]
        u_b = @view u[:, idx_b]
        
        phys = p.boundary_map[idx_a].physics #this doesn't impact performance at all really
        F = numerical_flux(phys, u_a, u_b, conn.area, conn.normal, conn.distance) #THIS DEFINITELY DOES
        #making F = 0.0 makes a 1 million cell sim take 2.6 seconds so numerical flux is the bottlneck
        
        @views du[:, idx_a] .-= F
        @views du[:, idx_b] .+= F
    end

    # Boundary & Source
    n_cells = size(u, 2)
    for i in 1:n_cells
        vol = p.mesh_cells[i].volume
        phys = p.boundary_map[i].physics

        bc = p.boundary_map[i]
        
        # Apply Capacity Division ( dU/dt = NetFlux / Volume )
        cap = capacity(phys, vol)
        @views du[:, i] ./= cap
        
        # Simple Boundary Condition Handling (Dirichlet Pinning)
        if bc.type == :Dirichlet
            @views du[:, i] .= 0.0
        end
    end
end

# ==============================================================================
# 4. SIMULATION: SOD SHOCK TUBE
# ==============================================================================


grid_dimensions = (1000, 1, 1) 
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 0.1, 0.1))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

# Define Physics (Air)
air_physics = CompressibleEulerPhysics(1.4, 287.0)

# Sets
addcellset!(grid, "high_pressure", x -> x[1] <= 0.5)
addcellset!(grid, "low_pressure", x -> x[1] > 0.5)

high_p_idxs = Set(getcellset(grid, "high_pressure"))
low_p_idxs = Set(getcellset(grid, "low_pressure"))

# Initial States (Sod Shock Tube standard values)
# Left: Rho=1, P=1, V=0
# Right: Rho=0.125, P=0.1, V=0

function get_conserved_state(rho, u, v, w, p, gamma)
    E = p / (gamma - 1.0) + 0.5 * rho * (u^2 + v^2 + w^2)
    return [rho, rho*u, rho*v, rho*w, E]
end

U_Left  = get_conserved_state(1.0, 0.0, 0.0, 0.0, 1.0, 1.4)
U_Right = get_conserved_state(0.125, 0.0, 0.0, 0.0, 0.1, 1.4)

function fluid_bc_mapper(grid, cell_id)
    # We use Neumann everywhere for the flux logic, but initialization sets the step
    return FluidBC(:Neumann, zeros(5), air_physics)
end

println("Building FVM System ($(grid_dimensions[1]) x $(grid_dimensions[2]) x $(grid_dimensions[3])) ($(grid_dimensions[1] * grid_dimensions[2] * grid_dimensions[3]) cells)")
@time fvm_prob = build_fvm_problem(grid, fluid_bc_mapper)
#a 1 million cell grid takes 3.5 seconds which is decent

# Initialize U0 (Matrix: 5 vars x N cells)
n_cells = length(fvm_prob.mesh_cells)
n_vars = 5 
u0 = zeros(Float64, n_vars, n_cells)

for i in 1:n_cells
    if i in high_p_idxs
        u0[:, i] = U_Left
    else
        u0[:, i] = U_Right
    end
end

# Time span is short because shocks move fast!
tspan = (0.0, 1.5)

# Using Implicit Solver (TRBDF2) for stability with the Jacobian detection
# Note: For pure shocks, Explicit methods (Tsit5) are often sharper, but TRBDF2 is robust.
detector = SparseConnectivityTracer.TracerSparsityDetector()
du0 = copy(u0)

#jac_sparsity = ADTypes.jacobian_sparsity(
    #(du, u) -> FVM_iter_f!(du, u, fvm_prob, 0.0), du0, u0, detector)

ode_func = ODEFunction(FVM_iter_f!)#, jac_prototype = float.(jac_sparsity))

prob = ODEProblem(ode_func, u0, tspan, fvm_prob)

println("Solving Navier-Stokes")
desired_amount_of_u = 100
@time sol = solve(prob, SSPRK43(), save_everystep=false)#, saveat=(tspan[end]/desired_amount_of_u), reltol=1e-3, abstol=1e-4) 
#for some reason, if I rerun sol = solve(...) with ctrl + enter to rerun it in the repl, I get a 1.417466 seconds (76 allocations: 1019.448 KiB) run time 
#WOAH, save_everystep=false makes a 4.1s 1.6M allocation go to a 1.5s 268k allocations for a 2000 cell problem

#we could use KYKSSPRK42() for DG in the future 
#1x1x1 takes 0.843465 seconds (1.58 M allocations: 85.729 MiB, 99.99% compilation time: 100% of which was recompilation)
#1000x1x1 takes 1.312537 seconds (1.60 M allocations: 87.932 MiB, 76.39% compilation time: 100% of which was recompilation)
#2000x1x1 takes 2.089783 seconds (1.60 M allocations: 88.545 MiB, 38.77% compilation time: 100% of which was recompilation)
#5000x1x1 takes 8.707805 seconds (1.60 M allocations: 89.965 MiB, 0.49% gc time, 10.30% compilation time: 100% of which was recompilation)
#strange that they all take around the same amount of allocations 

#that's probably good news because it means that there's probably a single source contributing to the ~1.6 M allocations
#based on the 1x1x1 example it seems like compilation takes a static 0.84 seconds

#Making the solution avaliable to paraview
record_sol = false

if record_sol
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

    for (step, t) in enumerate(sol.t)
        # We need to unpack the Matrix variables into separate vectors for Paraview
        
        # Get matrix at this timestep
        U_mat = sol.u[step] # Size (5, N)
        
        # Extract rows
        rho_field = vec(U_mat[1, :])
        E_field   = vec(U_mat[5, :])
        
        # Calculate Velocity Magnitude and Pressure for visualization
        vel_mag_field = zeros(n_cells)
        pressure_field = zeros(n_cells)
        
        for i in 1:n_cells
            col = @view U_mat[:, i]
            rho = col[1]
            mx, my, mz = col[2], col[3], col[4]
            E = col[5]
            
            # Velocity Mag
            v_mag = sqrt(mx^2 + my^2 + mz^2) / rho
            vel_mag_field[i] = v_mag
            
            # Pressure
            # p = (gamma-1)*(E - 0.5*rho*v^2)
            p = (1.4 - 1.0) * (E - 0.5 * rho * v_mag^2)
            pressure_field[i] = p
        end

        VTKGridFile(step_filename * "_$step.vtu", grid) do vtk
            write_cell_data(vtk, rho_field, "Density")
            write_cell_data(vtk, pressure_field, "Pressure")
            write_cell_data(vtk, vel_mag_field, "Velocity_Mag")
            pvd[t] = vtk
        end
    end
    vtk_save(pvd)
end