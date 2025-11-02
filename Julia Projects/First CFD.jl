using LinearAlgebra, SparseArrays, Plots, Printf, Unitful


rod_length = 10.0u"m"
T_left = 100.0u"K" 
T_right = 20.0u"K" 

rod_thermal_cond = 400u"W/(m*K)"
fluid_density = 10u"g*cm^-3"
fluid_specific_heat_capacity = 4.000u"J*g^-1*K^-1"
fluid_velocity = 3u"m/s"


function create_grid(L::Any, N::Int)
    unitful_dx = L / N
    dx = ustrip(unitful_dx)
    cell_centers = []
    face_coords = []
    for x in 1:N
        push!(cell_centers, (x - 0.5) * dx)
    end
    for x in 0:N
        push!(face_coords, x * dx)
    end
    return dx, cell_centers, face_coords
end

dx, cell_centers, face_coords = create_grid(rod_length, 10)

#remember dx is the length of each cell and k is the thermal conductivity 


function heat_flux_through_1d_grid(dx, cell_centers, face_coords, k, ρ, cp, velocity, T_left, T_right)
    N = length(cell_centers)
    T = Vector{Float64}(undef, N)

    dx = ustrip(dx)
    k = ustrip(k)
    ρ = ustrip(ρ)
    cp = ustrip(cp)
    velocity = ustrip(velocity)
    T_left = ustrip(T_left)
    T_right = ustrip(T_right)

    D = k / dx
    F = ρ * cp * velocity
    
    A = spzeros(N, N)
    b = zeros(N)

    for i in 1:N
        A[i, i] = 2.0 * D + abs(F)
        if i > 1 # Only if not the first cell
            A[i, i-1] = -(D + max(F, 0))
        end
        # East neighbor (T_{i+1})
        if i < N # Only if not the last cell
            A[i, i+1] = -(D + min(F, 0))
        end
        if i == 1 # First cell, West boundary
            # The T_A term (D * T_A) moves to the RHS of the equation
            b[i] += (D + max(F, 0))* T_left
        end
        if i == N # Last cell, East boundary
            # The T_B term (D * T_B) moves to the RHS of the equation
            b[i] += (D + min(F, 0)) * T_right
        end
    end
    return A, b
end

A, b = heat_flux_through_1d_grid(dx, cell_centers, face_coords, rod_thermal_cond, fluid_density, fluid_specific_heat_capacity, fluid_velocity, T_left, T_right)


println("\nMatrix A:")
display(A) # Use display for better sparse matrix output
println("\nVector b:")
display(b)

# 4. Solve the linear system
# Julia's backslash operator is highly optimized for this
T_solution_vals = A \ b
# Re-attach the unit to the solution vector
T_solution = T_solution_vals * Unitful.unit(T_left)

println("\nTemperature Solution (values): ", T_solution_vals)
println("Temperature Solution (with units): ", T_solution)

# 5. Verify and Visualize
# Analytical Solution: T(x) = T_A + (T_B - T_A) * x / L
# Make sure to work with stripped values for the analytical function for consistency
analytical_solution(x, L, T_A, T_B) = ustrip(T_A) + (ustrip(T_B) - ustrip(T_A)) * ustrip(x) / ustrip(L)

# Convert cell_centers to values for plotting and analytical solution
cell_centers_vals = ustrip.(cell_centers)
analytical_T = [analytical_solution(x, rod_length, T_left, T_right) for x in cell_centers_vals]


plot(cell_centers_vals, ustrip.(T_solution), label="FVM Solution", marker=:circle, line=:solid, linewidth=2,
     xlabel="Position (m)", ylabel="Temperature (K)", title="1D Steady-State Diffusion")
plot!(cell_centers_vals, analytical_T, label="Analytical Solution", linestyle=:dash, linewidth=2)
savefig("1d_diffusion_solution.png")
println("\nPlot saved to 1d_diffusion_solution.png")

#remember, you cannot evaluate this step by step, you're creating a relationship via 
#a matrix that is solved all at once 