using Ferrite
using LinearAlgebra
using SparseArrays
using NonlinearSolve
using RecursiveArrayTools
import SparseConnectivityTracer, ADTypes 
using SparseMatrixColorings
using AlgebraicMultigrid

# ==========================================
# 1. Physical Constants & Parameters
# ==========================================

const R_GAS = 8.314 # J/(mol K)

# 0: Solid, 1: Reforming Fluid, 2: Combustion Fluid
@enum ZoneType SOLID=0 REFORMING=1 COMBUSTION=2

struct PhysicsParams
    # -- Fluid Properties --
    mu_mix::Float64   # Viscosity [Pa.s]
    k_gas::Float64    # Thermal Cond Gas [W/mK]
    k_solid::Float64  # Thermal Cond Solid [W/mK]
    cp_gas::Float64   # Heat Capacity Gas [J/kgK]
    cp_solid::Float64 # Heat Capacity Solid [J/kgK]
    rho_solid::Float64
    D_mix::Float64    # Diffusivity [m^2/s]
    
    # -- Reaction Parameters (Xu-Froment / Simplified) --
    # Pre-exponential factors and Activation energies would go here
    # Using simplified placeholders for stability in this demo
    H_rxn_smr::Float64 # Enthalpy of SMR [J/mol] (Endothermic +)
    H_rxn_comb::Float64 # Enthalpy of Combustion [J/mol] (Exothermic -)
end

# Default values
params = PhysicsParams(
    3.0e-5, 0.05, 50.0, 2000.0, 500.0, 8000.0, 1.0e-5,
    206000.0, -802000.0 
)

# Variable Indices in the State Vector
const IDX_P = 1
const IDX_U = 2
const IDX_V = 3
const IDX_W = 4
const IDX_T = 5
const IDX_Y_CH4 = 6
const IDX_Y_H2O = 7
const IDX_Y_H2  = 8
const IDX_Y_CO  = 9
const IDX_Y_O2  = 10
# Y_CO2 is balance: 1 - sum(others)
const N_VARS = 10

# Molecular Weights (kg/mol)
const MW = Dict(
    :CH4 => 0.016, :H2O => 0.018, :H2 => 0.002, 
    :CO => 0.028, :O2 => 0.032, :CO2 => 0.044
)

# ==========================================
# 2. Data Structures
# ==========================================

struct CellData
    volume::Float64
    centroid::Vec{3, Float64}
    type::ZoneType
end

struct Connection
    cell_idx_a::Int
    cell_idx_b::Int
    area::Float64
    normal::Vec{3, Float64} # Normal points A -> B
    distance::Float64
    is_wall::Bool  # True if connection crosses Fluid/Solid boundary (Catalyst)
end

struct BoundaryFace
    cell_idx::Int
    area::Float64
    normal::Vec{3, Float64}
    type::Symbol # :Inlet, :Outlet, :Wall
end

struct FVMGeometry
    mesh_cells::Vector{CellData}
    connections::Vector{Connection}
    boundary_faces::Vector{BoundaryFace}
    # Indicies for easy looping
    fluid_idxs::Vector{Int}
    solid_idxs::Vector{Int}
end

# ==========================================
# 3. Grid Generation (Sandwich Structure)
# ==========================================

function build_geometry(dims::Tuple{Int,Int,Int})
    # Create a long channel: 3 layers in Y direction (Combustion, Steel, Reforming)
    # X: Flow direction
    L_x, L_y, L_z = 0.1, 0.03, 0.01
    
    # Standard Grid
    left = Vec{3}((0.0, 0.0, 0.0))
    right = Vec{3}((L_x, L_y, L_z))
    grid = generate_grid(Hexahedron, dims, left, right)
    
    n_cells = getncells(grid)
    cells_data = Vector{CellData}(undef, n_cells)
    
    # Y-thresholds for layers
    y_comb_end = L_y * 0.33
    y_ref_start = L_y * 0.66
    
    fluid_idxs = Int[]
    solid_idxs = Int[]

    # Geometry extraction helpers
    ref_shape = getrefshape(grid.cells[1])
    cell_qr = QuadratureRule{ref_shape}(2)
    poly_interp = Lagrange{ref_shape, 1}()
    cell_values = CellValues(cell_qr, poly_interp)
    
    for cell in CellIterator(grid)
        id = cellid(cell)
        Ferrite.reinit!(cell_values, cell)
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values))
        coords = getcoordinates(cell)
        cent = sum(coords) / length(coords)
        
        # Assign Types based on Y position
        y = cent[2]
        if y < y_comb_end
            z_type = COMBUSTION
            push!(fluid_idxs, id)
        elseif y > y_ref_start
            z_type = REFORMING
            push!(fluid_idxs, id)
        else
            z_type = SOLID
            push!(solid_idxs, id)
        end
        
        cells_data[id] = CellData(vol, cent, z_type)
    end
    
    # Connectivity
    top = ExclusiveTopology(grid)
    connections = Vector{Connection}()
    boundary_faces = Vector{BoundaryFace}()
    
    facet_qr = FacetQuadratureRule{ref_shape}(2)
    facet_values = FacetValues(facet_qr, poly_interp)
    
    for i in 1:n_cells
        for face_idx in 1:nfacets(grid.cells[i])
            neighbor_info = top.face_face_neighbor[i, face_idx]
            
            # Geometry of the face
            coords = getcoordinates(grid, i)
            Ferrite.reinit!(facet_values, coords, face_idx)
            area = sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values))
            n_ref = getnormal(facet_values, 1)

            if !isempty(neighbor_info)
                neighbor_idx = collect(neighbor_info[1].idx)[1]
                # Only process connection once
                if i < neighbor_idx
                    type_a = cells_data[i].type
                    type_b = cells_data[neighbor_idx].type
                    
                    cent_a = cells_data[i].centroid
                    cent_b = cells_data[neighbor_idx].centroid
                    dist = norm(cent_b - cent_a)
                    
                    # Detect Catalyst Wall (Interface between Solid and Fluid)
                    is_wall = (type_a == SOLID && type_b != SOLID) || (type_a != SOLID && type_b == SOLID)
                    
                    push!(connections, Connection(i, neighbor_idx, area, n_ref, dist, is_wall))
                end
            else
                # Boundary Face
                cent = cells_data[i].centroid
                # Simple BC classification
                bc_type = :Wall
                if cent[1] < 1e-4; bc_type = :Inlet; end
                if cent[1] > (L_x - 1e-4); bc_type = :Outlet; end
                
                push!(boundary_faces, BoundaryFace(i, area, n_ref, bc_type))
            end
        end
    end
    
    return FVMGeometry(cells_data, connections, boundary_faces, fluid_idxs, solid_idxs)
end

# ==========================================
# 4. Physics Helper Functions
# ==========================================

# Ideal Gas Law
function get_density(P, T, Y, type::ZoneType)
    if type == SOLID; return 8000.0; end
    
    # Calculate Mix MW
    # M_mix = 1 / sum(Y_i / M_i)
    # Reconstructing full Y vector including CO2
    Y_sum = Y[1]+Y[2]+Y[3]+Y[4]+Y[5]
    Y_CO2 = 1.0 - Y_sum
    
    inv_M = Y[1]/MW[:CH4] + Y[2]/MW[:H2O] + Y[3]/MW[:H2] + Y[4]/MW[:CO] + Y[5]/MW[:O2] + Y_CO2/MW[:CO2]
    M_mix = 1.0 / inv_M
    
    return (P * M_mix) / (R_GAS * T)
end

# Xu-Froment Kinetics (Simplified for Demo)
# Returns reaction rates [mol/(s*m2)] for SMR, WGS, SMR2
function xu_froment_rates(T, P, Y)
    # Real implementation requires partial pressures and adsorption constants
    # This is a functional placeholder to represent the load
    k1 = 1e8 * exp(-20000/T) 
    r_smr = k1 * Y[1] * Y[2] * P^2 # CH4 * H2O
    r_wgs = k1 * 0.1 * Y[4] * Y[2] * P^2 # CO * H2O
    
    # Coking Rate (Output only, simplified)
    r_coke = 1e-5 * Y[1] # Proportional to methane
    
    return [r_smr, r_wgs, 0.0], r_coke
end

# Combustion Kinetics
function combustion_rate(T, P, Y)
    # CH4 + 2O2 -> CO2 + 2H2O
    k = 1e9 * exp(-15000/T)
    return k * Y[1] * Y[5] * P^2 # CH4 * O2
end

# Upwind Scheme
function get_phi_interface(phi_L, phi_R, mass_flow)
    return (mass_flow > 0) ? phi_L : phi_R
end

# ==========================================
# 5. Residual Function (Monolithic FVM)
# ==========================================

function fvm_residual!(du, u, p, geo::FVMGeometry, phys::PhysicsParams)
    # 1. Unpack Dimensions
    n_cells = length(geo.mesh_cells)
    du .= 0.0
    
    # Helper to get variable at cell i
    # Layout: u[ (var-1)*n_cells + i ]
    function get_vars(idx)
        P = u[idx]
        vel = Vec{3}((u[n_cells + idx], u[2*n_cells + idx], u[3*n_cells + idx]))
        T = u[4*n_cells + idx]
        # Species
        Y = @view u[5*n_cells + idx : n_cells : 9*n_cells + idx] # 5 species
        return P, vel, T, Y
    end

    # 2. Loop Internal Connections (Fluxes)
    for conn in geo.connections
        iL, iR = conn.cell_idx_a, conn.cell_idx_b
        
        type_L = geo.mesh_cells[iL].type
        type_R = geo.mesh_cells[iR].type
        
        # Unpack States
        P_L, v_L, T_L, Y_L = get_vars(iL)
        P_R, v_R, T_R, Y_R = get_vars(iR)
        
        # --- A. Fluid Dynamics (Navier Stokes) ---
        # Only if both are fluid. If Solid-Solid -> Conduction only. If Solid-Fluid -> Wall.
        
        is_fluid_conn = (type_L != SOLID) && (type_R != SOLID)
        
        # Interpolate variables to face
        f_L = 0.5; f_R = 0.5 # Central differencing for diffusion properties
        
        # Normal vector n points L -> R
        n = conn.normal
        d = conn.distance
        A = conn.area
        
        # 1. Momentum & Mass Fluxes
        mass_flow = 0.0
        
        if is_fluid_conn
            rho_L = get_density(P_L, T_L, Y_L, type_L)
            rho_R = get_density(P_R, T_R, Y_R, type_R)
            rho_f = 0.5 * (rho_L + rho_R)
            
            # Rhie-Chow Interpolation (Simplified: Add pressure smoothing)
            # U_face = avg(U) - D * (grad P_face - avg(grad P))
            # Here: classic smoothing term proportional to (P_R - P_L)
            v_avg = 0.5 * (v_L + v_R)
            U_cont = dot(v_avg, n) - 0.01 * (P_R - P_L) 
            
            mass_flow = rho_f * U_cont * A
            
            # Continuity Residual
            du[iL] += mass_flow
            du[iR] -= mass_flow
            
            # Momentum Advection (Upwind)
            if mass_flow > 0
                v_up = v_L
            else
                v_up = v_R
            end
            mom_flux = mass_flow * v_up
            
            # Pressure Force
            p_force = 0.5 * (P_L + P_R) * A * n
            
            # Viscous Stress (Simplified Laplacian)
            mu = phys.mu_mix
            tau = mu * (v_R - v_L) / d * A
            
            # Momentum Residuals (x, y, z)
            for k in 1:3
                idx_v = k + 1
                du[(idx_v-1)*n_cells + iL] += mom_flux[k] + p_force[k] - tau[k]
                du[(idx_v-1)*n_cells + iR] -= mom_flux[k] + p_force[k] - tau[k]
            end
        end
        
        # --- B. Heat Transfer ---
        # Thermal Conductivity
        k_L = (type_L == SOLID) ? phys.k_solid : phys.k_gas
        k_R = (type_R == SOLID) ? phys.k_solid : phys.k_gas
        k_eff = 2*k_L*k_R / (k_L + k_R)
        
        # Conduction
        q_cond = -k_eff * (T_R - T_L) / d * A
        
        # Advection (Enthalpy)
        cp = phys.cp_gas # Simplified constant Cp
        h_up = cp * ((mass_flow > 0) ? T_L : T_R)
        q_adv = mass_flow * h_up
        
        # Energy Residual
        du[4*n_cells + iL] += q_adv + q_cond
        du[4*n_cells + iR] -= q_adv + q_cond
        
        # --- C. Species Transport ---
        if is_fluid_conn
            D = phys.D_mix
            rho_f = 0.5 * (get_density(P_L, T_L, Y_L, type_L) + get_density(P_R, T_R, Y_R, type_R))
            
            for s in 1:5 # Loop 5 tracked species
                # Diffusion
                J_diff = -rho_f * D * (Y_R[s] - Y_L[s]) / d * A
                
                # Advection
                Y_up = (mass_flow > 0) ? Y_L[s] : Y_R[s]
                J_adv = mass_flow * Y_up
                
                idx_s = 4 + s
                du[idx_s*n_cells + iL] += J_adv + J_diff
                du[idx_s*n_cells + iR] -= J_adv + J_diff
            end
        end
        
        # --- D. Wall / Catalyst Reaction ---
        if conn.is_wall
            # Identify which is fluid and which is solid
            idx_f, idx_s = (type_L != SOLID) ? (iL, iR) : (iR, iL)
            is_reforming = (geo.mesh_cells[idx_f].type == REFORMING)
            
            P_f, _, T_f, Y_f = get_vars(idx_f)
            
            if is_reforming
                # SMR Kinetics (Xu Froment)
                rates, _ = xu_froment_rates(T_f, P_f, Y_f)
                r1, r2, r3 = rates[1], rates[2], rates[3]
                
                # Stoichiometry (mol) * Rate (mol/s) * MW (kg/mol)
                # SMR: CH4 + H2O -> CO + 3H2
                # WGS: CO + H2O -> CO2 + H2
                
                S_CH4 = (-r1 - r3) * MW[:CH4] * A
                S_H2O = (-r1 - r2 - 2*r3) * MW[:H2O] * A
                S_H2  = (3*r1 + r2 + 4*r3) * MW[:H2] * A
                S_CO  = (r1 - r2) * MW[:CO] * A
                # CO2 is balance, not solved explicitly
                
                # Apply Source to Fluid Cell
                du[5*n_cells + idx_f] -= S_CH4
                du[6*n_cells + idx_f] -= S_H2O
                du[7*n_cells + idx_f] -= S_H2
                du[8*n_cells + idx_f] -= S_CO
                
                # Heat of Reaction (Applied to Solid Interface or Split?)
                # SMR is Endothermic (consumes heat). Source is Negative.
                # Q = Rate * DeltaH
                Q_rxn = -(r1 * phys.H_rxn_smr) * A
                
                # Apply heat sink to Fluid (or split with solid)
                du[4*n_cells + idx_f] -= Q_rxn 
            end
        end
    end
    
    # 3. Volumetric Reactions (Combustion) & BCs
    for i in 1:n_cells
        cell = geo.mesh_cells[i]
        vol = cell.volume
        
        # Volumetric Combustion
        if cell.type == COMBUSTION
            P, _, T, Y = get_vars(i)
            r_comb = combustion_rate(T, P, Y) # mol/m3 s
            
            # CH4 + 2O2 -> CO2 + 2H2O
            S_CH4 = -r_comb * MW[:CH4] * vol
            S_O2  = -2 * r_comb * MW[:O2] * vol
            S_H2O = 2 * r_comb * MW[:H2O] * vol
            
            du[5*n_cells + i] -= S_CH4
            du[9*n_cells + i] -= S_O2
            du[6*n_cells + i] -= S_H2O
            
            # Exothermic Heat (Source Positive)
            Q_comb = r_comb * abs(phys.H_rxn_comb) * vol
            du[4*n_cells + i] -= Q_comb # Residual = Flux_out - Source -> so minus
        end
        
        # Solid Domain Constraints (Fix Velocity = 0)
        if cell.type == SOLID
            du[i] = 0.0 # Pressure (undefined, set to 0 or 1 atm ref)
            du[n_cells + i] = u[n_cells + i]     # u = 0
            du[2*n_cells + i] = u[2*n_cells + i] # v = 0
            du[3*n_cells + i] = u[3*n_cells + i] # w = 0
            # Species = 0
            for s in 5:9
                du[s*n_cells + i] = u[s*n_cells + i]
            end
        end
    end
    
    # 4. Boundary Conditions
    for bc in geo.boundary_faces
        i = bc.cell_idx
        # Simple Dirichlet implementation via penalty or direct replacement
        # Here we add flux contribution for Inlet/Outlet
        
        if bc.type == :Inlet
            # Fixed Inlet State
            if geo.mesh_cells[i].type == REFORMING
                # SMR Inlet: 500C, 2 bar, High CH4/H2O
                T_in = 773.0; P_in = 2e5; vel_in = Vec{3}((1.0, 0.0, 0.0))
                Y_in = [0.2, 0.6, 0.0, 0.0, 0.0] # CH4, H2O
            elseif geo.mesh_cells[i].type == COMBUSTION
                # Comb Inlet: 1000C, 1 bar, CH4/O2
                T_in = 1273.0; P_in = 1e5; vel_in = Vec{3}((2.0, 0.0, 0.0))
                Y_in = [0.1, 0.0, 0.0, 0.0, 0.2] # CH4, ..., O2
            else
                continue # Solid inlet??
            end
            
            # Add fluxes based on fixed inlet values
            # (Simplified: forcing residual to match value)
            # A robust code would compute flux (rho u A)_in * (phi_cell - phi_in)
            
            # Force Dirichlet (Soft)
            penalty = 1e5
            du[4*n_cells + i] += penalty * (u[4*n_cells + i] - T_in)
            du[i] += penalty * (u[i] - P_in)
            # ... and so on for species
        elseif bc.type == :Outlet
            # Fixed Pressure
            P_out = 101325.0
            du[i] += 1e5 * (u[i] - P_out)
        elseif bc.type == :Wall
            # No Slip (Velocity = 0)
            for k in 1:3
                du[k*n_cells + i] += 1e6 * u[k*n_cells + i]
            end
            # Adiabatic walls (Flux=0) is implicitly handled by doing nothing
        end
    end
end

# ==========================================
# 6. Main Execution Setup
# ==========================================

println("Generating SMR Microchannel Geometry...")
# 5x10x2 grid for demonstration (Sandwich structure needs Y resolution)
# 10 cells in X (flow), 9 in Y (3 Comb, 3 Steel, 3 Ref), 2 in Z
mesh_dims = (10, 9, 2) 
geo = build_geometry(mesh_dims)
n_cells = length(geo.mesh_cells)

println("Initializing State Vector ($n_cells cells, $N_VARS vars/cell)...")
# Flattened state vector u
# Order: [P_1...P_n, u_1...u_n, v_1...v_n, w_1...w_n, T_1...T_n, Y1_1...Y1_n, ...]
u0 = zeros(n_cells * N_VARS)

for i in 1:n_cells
    type = geo.mesh_cells[i].type
    
    # Defaults
    P = 101325.0
    T = 800.0
    vel = 0.1
    Y = [0.0, 0.0, 0.0, 0.0, 0.0]
    
    if type == REFORMING
        P = 2e5
        T = 773.0
        vel = 1.0
        Y = [0.2, 0.6, 0.01, 0.0, 0.0] # CH4, H2O
    elseif type == COMBUSTION
        T = 1200.0
        vel = 2.0
        Y = [0.1, 0.05, 0.0, 0.0, 0.2] # CH4, O2
    else # Solid
        T = 900.0
        vel = 0.0
    end
    
    u0[i] = P
    u0[n_cells + i] = vel # u_x
    u0[4*n_cells + i] = T
    
    for s in 1:5
        u0[(4+s)*n_cells + i] = Y[s]
    end
end

# Define Problem
f_closure = (du, u, p) -> fvm_residual!(du, u, p, geo, params)

# Sparse Connectivity Detection
println("Detecting Sparsity Pattern...")
detector = SparseConnectivityTracer.TracerLocalSparsityDetector()
du0 = similar(u0) .*= 0.0
jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> f_closure(du, u, nothing), du0, u0, detector)

println("Sparsity density: ", nnz(jac_sparsity) / length(jac_sparsity))


coloring_prob = ColoringProblem(partition = :column)
coloring_alg = GreedyColoringAlgorithm(; decompression = :direct)
coloring_result = coloring(jac_sparsity, coloring_prob, coloring_alg)
colors = column_colors(coloring_result)

nl_func = NonlinearFunction(f_closure, jac_prototype = float.(jac_sparsity))

prob = NonlinearProblem(nl_func, u0)

println("Solving SMR System (Newton-Raphson)...")
# Using a robust linear solver (GMRES) inside Newton because the system is non-symmetric and stiff
# relax damping to help convergence


function algebraicmultigrid(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(convert(AbstractMatrix, W)))
    else
        Pl = Plprev
    end
    Pl, nothing
end

sol = solve(prob, NewtonRaphson(linsolve = KrylovJL_GMRES()))
println(sol.u)
#linesearch=LineSearchesJL()