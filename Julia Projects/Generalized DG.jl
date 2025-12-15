using DifferentialEquations
using LinearAlgebra #for norm()
using WriteVTK
using Dates
#using Unitful
using SparseArrays
using RecursiveArrayTools
import SparseConnectivityTracer, ADTypes 
using ProfileView
using SparseMatrixColorings
using ILUZero
using Ferrite
import AlgebraicMultigrid

abstract type AbstractPhysics end

struct CellData
    mass_matrix::AbstractMatrix #to be refined later, don't know what kind of matrix I'll use right now
    stiffness_matrix::AbstractMatrix
    centroid::Vec{3, Float64}
    N_points::Int
    polynomial_func::Function
    polynomial_order
end

struct Connection
    cell_idx_a::Int
    cell_idx_b::Int
    area::Float64
    normal::Vec{3, Float64}
    face_quadrature_points::Vector{Float64}
end

abstract type AbstractBC end

struct FVMProblem{B<:AbstractBC}
    mesh_cells::Vector{CellData}
    connections::Vector{Connection}
    boundary_map::Vector{B}
    free_idxs::Vector{Int} 
    dirichlet_idxs::Vector{Int}
end

struct HeatPhysics <: AbstractPhysics
    k::Float64 
    rho::Float64
    cp::Float64
    source_term::Float64 #volumetric heating
end

struct HeatBC <: AbstractBC
    type::Symbol #dirichlet or neumann
    initial::Float64
    physics::HeatPhysics
    N_points::Int
    polynomial_func::Function
end

function numerical_flux(phys::HeatPhysics, T_L, T_R, area, normal, dist)
    grad_T = (T_R - T_L) / dist

    q = -phys.k * grad_T

    return q * area
end

function source(phys::HeatPhysics, u, vol)
    return phys.source_term * vol
end

function capacity(phys::HeatPhysics, vol)
    return phys.rho * phys.cp * vol
end


function build_fvm_problem(grid, bc_map_func)
    n_cells = getncells(grid)
    
    ref_shape = getrefshape(grid.cells[1])
    poly_interp = Lagrange{ref_shape, 1}() #1 = linear elements 2 = quadratic/curved edges
    facet_qr = FacetQuadratureRule{ref_shape}(2) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation
    facet_values = FacetValues(facet_qr, poly_interp)
    cell_qr = QuadratureRule{ref_shape}(2)
    cell_values = CellValues(cell_qr, poly_interp)

    cells_data = Vector{CellData}(undef, n_cells)
    
    #volumes and centroids
    for cell in CellIterator(grid)
        id = cellid(cell)
        Ferrite.reinit!(cell_values, cell)
        vol = sum(getdetJdV(cell_values, qp) for qp in 1:getnquadpoints(cell_values))
        coords = getcoordinates(cell)
        cent = sum(coords) / length(coords)
        boundary_map[i] = bc_map_func(grid, i)
        cells_data[id] = CellData(vol, cent)
    end

    top = ExclusiveTopology(grid)
    connections = Vector{Connection}()
    
    #store bcs
    boundary_map = Vector{HeatBC}(undef, n_cells)

    free_idxs = Vector{Int}()
    dirichlet_idxs = Vector{Int}()
    sizehint!(free_idxs, n_cells)
    sizehint!(dirichlet_idxs, floor(Int, n_cells/10))
    
    for i in 1:n_cells
        boundary_map[i] = bc_map_func(grid, i)

        Ferrite.reinit!(cell_values, cell)

        if boundary_map[i].type == :Dirichlet
            push!(dirichlet_idxs, i)
        else
            push!(free_idxs, i)
        end

        for face_idx in 1:nfacets(grid.cells[i])
            neighbor_info = top.face_face_neighbor[i, face_idx]
            
            if !isempty(neighbor_info)
                neighbor_idx = collect(neighbor_info[1].idx)[1]
                
                #avoid duplicates
                if i < neighbor_idx
                    coords = getcoordinates(grid, i)
                    Ferrite.reinit!(facet_values, coords, face_idx)
                    
                    area = sum(getdetJdV(facet_values, qp) for qp in 1:getnquadpoints(facet_values))
                    
                    n_ref = getnormal(facet_values, 1) 
                    
                    #centroid_a = cells_data[i].centroid
                    #centroid_b = cells_data[neighbor_idx].centroid
                    #dist = norm(centroid_b - centroid_a)

                    push!(connections, Connection(i, neighbor_idx, area, n_ref, getweights(grid.cells[i], face_idx)))
                end
            end
        end
    end
    
    return FVMProblem(cells_data, connections, boundary_map, free_idxs, dirichlet_idxs)
end

function FVM_iter_f!(du, u, p::FVMProblem, t)
    du .= 0.0

    #internal flux loop
    for conn in p.connections
        idx_a = conn.cell_idx_a
        idx_b = conn.cell_idx_b
        
        #perhaps it would be quicker to evaluate this polynomial at the beginning to prevent b's polynomial from being calculated twice if it had two connections
        u_a_vec = @view u[idx_a, :]
        a_val = p.boundary_map[i].polynomial_func(u_a_vec)
        u_b_vec = @view u[idx_b, :]
        b_val = p.boundary_map[idx_b].polynomial_func(u_b_vec)
        
        F = numerical_flux(p.boundary_map[idx_a].physics, a_val, b_val, conn.area, conn.normal, conn.distance)

        @views du[idx_a, :] .-= F
        @views du[idx_b, :] .+= F
    end

    #bcs loop
    for i in p.free_idxs
        vol = p.mesh_cells[i].volume

        u_curr = @view u[i, :]
        S = source(p.boundary_map[i].physics, u_curr, vol)

        @views du[i, :] .+= S
        
        cap = capacity(p.boundary_map[i].physics, vol)
        @views du[i, :] ./= cap
    end

    for i in p.dirichlet_idxs
        du[1, i] = 0.0
    end
end


grid_dimensions = (200, 10, 10)
#grid_dimensions = (100, 100, 100)
left = Ferrite.Vec{3}((0.0, 0.0, 0.0))
right = Ferrite.Vec{3}((1.0, 1.0, 1.0))
grid = generate_grid(Hexahedron, grid_dimensions, left, right)

cell_half_dist = (right[1] / grid_dimensions[1]) / 2 #this is to deal with uneven x grid dimensions
addcellset!(grid, "copper", x -> x[1] <= (left[1] + (right[1] / 2) + cell_half_dist))
getcellset(grid, "copper")
addcellset!(grid, "steel", x -> x[1] >= left[1] + (right[1] / 2))
getcellset(grid, "steel")

heated_copper_physics = HeatPhysics(401.0, 8960.0, 385.0, 0)#80000.0) #note that if we created a new struct for copper_physics performance would die
steel_physics = HeatPhysics(30.0, 8000.0, 460.0, 0)

global copper_cell_set_idxs = Set(getcellset(grid, "copper"))
global steel_cell_set_idxs = Set(getcellset(grid, "steel"))

function poly_func(u_vec)
    return u_vec[1]
end

function my_bc_mapper(grid, cell_id)
    if cell_id in copper_cell_set_idxs
        return HeatBC(:Neumann, 500.0, heated_copper_physics, 8, poly_func) #use :Dirichlet to fix temperature to initial in HeatBC
    elseif cell_id in steel_cell_set_idxs
        return HeatBC(:Neumann, 300.0, steel_physics, 8, poly_func)
    end
end

println("fvm system build time")
#@time fvm_prob = build_fvm_problem(grid, my_bc_mapper)

cell_weights_vec = Vector[]
face_weights_vec = Vector[]

ref_shape = getrefshape(grid.cells[1])
poly_interp = Lagrange{ref_shape, 1}() #1 = linear elements 2 = quadratic/curved edges
facet_qr = FacetQuadratureRule{ref_shape}(1) #2 represents the number of integration points. Basically higher number = higher accuracy but more computation
facet_values = FacetValues(facet_qr, poly_interp)
cell_qr = QuadratureRule{ref_shape}(2)
cell_values = CellValues(cell_qr, poly_interp)

for cell in CellIterator(grid)
    i = cellid(cell)
    Ferrite.reinit!(cell_values, cell)
    push!(cell_weights_vec, Ferrite.getweights(cell_qr))
    face_weights_sub_vec = Vector[]
    for face_idx in 1:nfacets(grid.cells[i])
        push!(face_weights_sub_vec, Ferrite.getweights(facet_qr, face_idx))
    end
    push!(face_weights_vec, face_weights_sub_vec)
end

cell_weights_vec
face_weights_vec



n_cells = length(fvm_prob.mesh_cells)
n_vars = 1
n_points = 8

u0 = zeros(Float64, n_cells, n_vars, n_points)

for i in 1:n_cells
    u0[i, 1, :] .= fvm_prob.boundary_map[i].initial
end

# THIS IS FOR IMPLICIT
#
tspan = (0.0, 100000.0)

detector = SparseConnectivityTracer.TracerSparsityDetector()
du0 = copy(u0)

jac_sparsity = ADTypes.jacobian_sparsity(
    (du, u) -> FVM_iter_f!(du, u, fvm_prob, 0.0), du0, u0, detector)

coloring_prob = ColoringProblem(partition = :column)
coloring_alg = GreedyColoringAlgorithm(; decompression = :direct)
coloring_result = coloring(jac_sparsity, coloring_prob, coloring_alg)
colors = column_colors(coloring_result)
ncolors(coloring_result)

#comp_jac_sparsity = compress(jac_sparsity, colors)

ode_func = ODEFunction(FVM_iter_f!, jac_prototype = float.(jac_sparsity), colorvec=colors)
#jac_prototype = float.(jac_sparsity))#, colorvec = colors.color)#, jac_prototype = float.(jac_sparsity))
#adding the colorvec = colors.color doesn't do anything
#ode_func = ODEFunction(FVM_iter_f!, colorvec = colors.color) 
prob = ODEProblem(ode_func, u0, tspan, fvm_prob)

#future notes
    # for sparsity detection use SparseConnectivityTracer
    # for coloring use SparseMatrixColorings
    # for differentiation use DifferentiationInterface.jl
    # see here for more: https://discourse.julialang.org/t/sparsedifftools-polyesterforwarddiff-pre-allocation-difficulties/130292/6

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

desired_amount_of_u = 100

#@time sol = solve(prob, FBDF(linsolve = KrylovJL_GMRES(), precs = iluzero, concrete_jac = true), saveat=(tspan[end]/desired_amount_of_u)) #save_everystep = false)
#also check out HYPRE.jl as it has PCG

#@time sol = solve(prob, FBDF(linsolve = KLUFactorization(reuse_symbolic=true)), saveat=(tspan[end]/desired_amount_of_u)) #save_everystep = false)

#@time sol = solve(prob, FBDF(linsolve = KrylovJL_GMRES(precs = WeightedDiagonalPreconBuilder(w = 0.9))), saveat=(tspan[end]/desired_amount_of_u)) #save_everystep = false)

# VERY IMPORTANT: I just learned that Ctrl+Enter is a better representation of speed than rerunning the entire file with ctrl+n because structs are re-initialized
sol.destats
#KLUFactorization(reuse_symbolic=true)) barely changes anything
#KrylovJL_GMRES() is very, very slow

#@time sol = solve(prob, FBDF(linsolve = KrylovJL_GMRES(precs=ADTypes.NoPreconditioner(), concrete_jac=true)), saveat=(tspan[end]/desired_amount_of_u)) #save_everystep = false)
#we'll probably want to switch to the above later
#ILUZero is a cheap preconditioner

#prob = ODEProblem(FVM_iter_f!, u0, tspan, fvm_prob)
#

#THIS IS FOR EXPLICIT
#=
tspan = (0.0, 100.0)
prob = ODEProblem(FVM_iter_f!, u0, tspan, fvm_prob)

println("sol time")
desired_amount_of_u = 100
@time solve(prob, Tsit5(), tspan=(0.0, 0.001), saveat=(tspan[end]/1)) 
@time sol = solve(prob, Tsit5(), saveat=(tspan[end]/desired_amount_of_u)) 
#whaaah, for some reason whenever I change to Float32 instead of Float64 it is shorter but it takes more memory and more allocations? What?!?
#Float32 run: 0.755455 seconds (702 allocations: 36.193 MiB, 2.22% gc time)
#Float64 run: 0.100255 seconds (405 allocations: 27.028 MiB) #where is the garbage collection coming from?
#Also, whether the garbage collector is triggered or not seems very inconsistent 
#on Float64 it happens sometimes
#on Float32 it happens everytime
    #sometimes on Float32 GC can take up 50% of total time!
=#
sol.destats

#VSCodeServer.@profview
#@time
#saveat=(tspan[end]/desired_amount_of_u)) 
#save_everystep = false

#Making the solution avaliable to paraview
record_sol = false

length(sol.u[1, :])

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
        temp_field = sol.u[step][:, 1]
        VTKGridFile(step_filename * " $step" * " at $(t)s.vtu", grid) do vtk
            write_cell_data(vtk, temp_field, "T")
            pvd[t] = vtk
        end
    end
    vtk_save(pvd)
end
