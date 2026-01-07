using Unitful
using CoolProp
using Roots

function reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, characteristic_length)
    fluid_density = PropsSI("D", "T", fluid_temp, "P", fluid_pressure, fluid)
    dynamic_viscosity = PropsSI("V", "T", fluid_temp, "P", fluid_pressure, fluid)
    return (fluid_density * fluid_velocity * characteristic_length) / (dynamic_viscosity)
end
    
function prandtl_number(fluid, fluid_temperature, fluid_pressure)
    dynamic_viscosity = PropsSI("V", "T", fluid_temperature, "P", fluid_pressure, fluid)
    thermal_cond = PropsSI("conductivity", "T", fluid_temperature, "P", fluid_pressure, fluid)
    specific_heat = PropsSI("C", "T", fluid_temperature, "P", fluid_pressure, fluid)
    density = PropsSI("D", "T", fluid_temperature, "P", fluid_pressure, fluid)
    return (dynamic_viscosity * specific_heat) / (thermal_cond)
end

prandtl_number("water", 300u"K", 1u"bar")

function grashof_number(ambient_fluid, ambient_fluid_temperature, ambient_fluid_pressure, characteristic_length, heat_exchange_surface_temperature)
    gravitational_acceleration = 9.81u"m/s"
    coefficient_of_expansion = PropsSI("isobaric_expansion_coefficient", "T", ambient_fluid_temperature, "P", ambient_fluid_pressure, ambient_fluid)
    
    dynamic_viscosity = PropsSI("V", "T", ambient_fluid_temperature, "P", ambient_fluid_pressure, ambient_fluid)
    density = PropsSI("D", "T", ambient_fluid_temperature, "P", ambient_fluid_pressure, ambient_fluid)
    kinematic_viscosity = dynamic_viscosity / density
    
    return (gravitational_acceleration * coefficient_of_expansion * (heat_exchange_surface_temperature - ambient_fluid_temperature) * (characteristic_length^3)) / (kinematic_viscosity^2)
end

function laminar_entry_thermal_Hausen(Re, Pr, pipe_length, pipe_diameter)
    Gz = (pipe_diameter / pipe_length) * Re * Pr
    return 3.66 + (0.0668*Gz)/(1+0.04*(Gz)^(2/3.0))
end

function ht_internal_conv(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_diameter, pipe_length, absolute_pipe_roughness)
    Re = reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_diameter)
    Pr = prandtl_number(fluid, fluid_temp, fluid_pressure)
    return laminar_entry_thermal_Hausen(Re, Pr, pipe_length, pipe_diameter)
end

function Nu_cylinder_Sanitjai_Goldstein(Re, Pr)
    return  0.446 * Re^0.5 * Pr^0.35 + 0.528 * ((6.5^-5 * exp(-5 * Re/5000.0))
    + (0.031 * Re^0.8)^-5)^-0.2 * Pr^0.42
end

function ht_external_pipe(external_fluid, external_fluid_temp, external_fluid_pressure, fluid_velocity, pipe_diameter)
    Re = reynolds_number(external_fluid, external_fluid_temp, external_fluid_pressure, fluid_velocity, pipe_diameter)
    
    PR = prandtl_number(external_fluid, external_fluid_temp, external_fluid_pressure)

    return Nu_cylinder_Sanitjai_Goldstein(Re, PR)
end

function overall_heat_transfer_coeff(
    internal_nusselt_number, internal_fluid_cond, internal_characteristic_length, internal_area, 
    external_nusselt_number, external_fluid_cond, external_characteristic_length, external_area, 
    layers)

    #Keep in mind that layers if formatted [(thickness, thermal conductivity, area), (thickness, thermal conductivity, area)]
    conductive_resistance = sum([layer.thickness / (layer.thermal_conductivity * layer.area) for layer in layers])
    
    internal_heat_transfer_coeff = (internal_nusselt_number * internal_fluid_cond) / internal_characteristic_length
    external_heat_transfer_coeff = (external_nusselt_number * external_fluid_cond) / external_characteristic_length
    
    return 1 / ((1 / (internal_heat_transfer_coeff * internal_area)) + (conductive_resistance) + (1 / (external_heat_transfer_coeff * external_area)))
end

function get_overall_heat_transfer_coefficent(
    internal_fluid, internal_fluid_temp, internal_fluid_pressure, internal_fluid_velocity,
    external_fluid, external_fluid_temp, external_fluid_pressure, external_fluid_velocity,
    pipe_inner_diameter, pipe_outer_diameter, 
    pipe_length, internal_pipe_absolute_roughness,
    layers
    )
    
    internal_area = (pi * pipe_inner_diameter) * pipe_length
    internal_thermal_cond = PropsSI("conductivity", "T", internal_fluid_temp, "P", internal_fluid_pressure, internal_fluid)
    internal_NU = ht_internal_conv(
        internal_fluid, internal_fluid_temp, internal_fluid_pressure, internal_fluid_velocity, 
        pipe_inner_diameter, pipe_length, internal_pipe_absolute_roughness
    )
    
    external_area = (pi * pipe_inner_diameter) * pipe_length
    external_thermal_cond = PropsSI("conductivity", "T", external_fluid_temp, "P", external_fluid_pressure, external_fluid)
    external_NU = ht_external_pipe(
        external_fluid, external_fluid_temp, external_fluid_pressure, external_fluid_velocity, 
        pipe_outer_diameter
    ) 

    layers[1].area = (pi * pipe_inner_diameter) * pipe_length
    
    #Keep in mind that layers if formatted [(thickness, thermal conductivity, area), (thickness, thermal conductivity, area)]
    overall_heat_trans_coeff = overall_heat_transfer_coeff(
        internal_NU, internal_thermal_cond, pipe_inner_diameter, internal_area, 
        external_NU, external_thermal_cond, pipe_outer_diameter, external_area, 
        layers
    )
    return overall_heat_trans_coeff
end

mutable struct Layer
    thickness
    thermal_conductivity
    area
end

pipe_length = 1.0u"m"
pipe_inner_diameter = 10u"cm"
pipe_outer_diameter = 11u"cm"
thickness = pipe_outer_diameter - pipe_inner_diameter
pipe_area = (pi * pipe_inner_diameter) * pipe_length

copper_layer = Layer(thickness, 400u"W/(m*K)", pipe_area)

h_target = 0.5u"W/K"
solve_for_pipe_diameter = x -> (get_overall_heat_transfer_coefficent(
    "water", #internal fluid
    300u"K", #internal temp
    1u"bar", #internal pressure
    1u"m/s", #internal fluid velocity

    "water", #external fluid
    450u"K", #external temp
    1u"bar", #external pressure
    1u"m/s", #external fluid velocity

    x * 1.0u"cm", #pipe ID
    (x+1) * 1.0u"cm", #pipe OD
    1u"m", #pipe length
    0.0015u"mm", #pipe absolute roughness
    [copper_layer]
) - h_target)

find_zero(solve_for_pipe_diameter, 0.5)

function min_heat_cap_rate(cp_hot, m_hot, cp_cold, m_cold)
    heat_cap_rate_hot = cp_hot * m_hot  
    heat_cap_rate_cold = cp_cold * m_cold   
    return min(heat_cap_rate_hot, heat_cap_rate_cold)
end

function heat_cap_ratio(cp_hot, m_hot, cp_cold, m_cold)
    heat_cap_rate_hot = cp_hot * m_hot  
    heat_cap_rate_cold = cp_cold * m_cold   
    cp_min = min(heat_cap_rate_hot, heat_cap_rate_cold)
    cp_max = max(heat_cap_rate_hot, heat_cap_rate_cold)
    return cp_min/cp_max
end

function effectiveness_from_NTU(NTU, Cr)
    return (1.0 - exp(-NTU*(1.0 - Cr)))/(1.0 - Cr*exp(-NTU*(1.0 - Cr)))
end

function effectiveness_NTU_method(
    hot_inlet_temp, cp_hot, m_hot, 
    cold_inlet_temp, cp_cold, m_cold, 
    UA)
    
    cp_min = min_heat_cap_rate(cp_hot, m_hot, cp_cold, m_cold)

    NTU = UA / cp_min

    Cr = heat_cap_ratio(cp_hot, m_hot, cp_cold, m_cold)

    effectiveness = effectiveness_from_NTU(NTU, Cr)    

    Q = effectiveness * cp_min * (hot_inlet_temp - cold_inlet_temp)

    return Q
end

function total_heat_duty(
    hot_inlet_temp, m_hot, 
    cold_inlet_temp, m_cold,
    
    internal_fluid, internal_fluid_temp, internal_fluid_pressure, internal_fluid_velocity,
    external_fluid, external_fluid_temp, external_fluid_pressure, external_fluid_velocity,
    pipe_inner_diameter, pipe_outer_diameter, 
    pipe_length, internal_pipe_absolute_roughness,
    layers
    )

    overall_heat_trans_coeff = get_overall_heat_transfer_coefficent(
        internal_fluid, internal_fluid_temp, internal_fluid_pressure, internal_fluid_velocity,
        external_fluid, external_fluid_temp, external_fluid_pressure, external_fluid_velocity,
        pipe_inner_diameter, pipe_outer_diameter, 
        pipe_length, internal_pipe_absolute_roughness,
        layers
    )

    hot_fluid_specific_heat = PropsSI("C", "T", internal_fluid_temp, "P", internal_fluid_pressure, internal_fluid)
    cold_fluid_specific_heat = PropsSI("C", "T", external_fluid_temp, "P", external_fluid_pressure, external_fluid)

    Q = effectiveness_NTU_method(
        hot_inlet_temp, hot_fluid_specific_heat, m_hot, 
        cold_inlet_temp, cold_fluid_specific_heat, m_cold, 
        overall_heat_trans_coeff
    )
    
    return Q
end 

Q_per_m_target = 4000.0u"W/m"

hot_fluid_inlet_temp = 350.0u"K"
hot_fluid_m = 1.0u"kg/s"

cold_fluid_inlet_temp = 300.0u"K"
cold_fluid_m = 1.0u"kg/s"

find_pipe_diameter_from_heat_duty = pipe_diameter -> ((total_heat_duty(
    hot_fluid_inlet_temp, hot_fluid_m, 
    cold_fluid_inlet_temp, cold_fluid_m,
    "water", #internal fluid
    hot_fluid_inlet_temp, #internal tempm
    1.0u"bar", #internal pressure
    1.0u"m/s", #internal fluid velocity

    "water", #external fluid
    cold_fluid_inlet_temp, #external temp
    1.0u"bar", #external pressure
    1.0u"m/s", #external fluid velocity

    pipe_diameter * 1.0u"cm", #pipe ID
    (pipe_diameter + 1) * 1.0u"cm", #pipe OD
    1.0u"m", #pipe length
    0.0015u"mm", #pipe absolute roughness
    [copper_layer]
) / pipe_length) - Q_per_m_target)

find_zero(find_pipe_diameter_from_heat_duty, 0.5) * 1.0u"cm"
#returns 5.1233 cm which seems reasonable

desired_iters = 2

max_Q_per_m = 100000u"W/m"

interval = max_Q_per_m / desired_iters

Q_per_m_span = 200.0u"W/m":interval:max_Q_per_m

Q_per_m_vals = []

corresponding_pipe_diameter_vals = []
#=
for (i, Q) in enumerate(Q_per_m_span)
    find_pipe_diameter_from_heat_duty = pipe_diameter -> ((total_heat_duty(
        hot_fluid_inlet_temp, hot_fluid_m, 
        cold_fluid_inlet_temp, cold_fluid_m,
        "water", #internal fluid
        hot_fluid_inlet_temp, #internal tempm
        1.0u"bar", #internal pressure
        1.0u"m/s", #internal fluid velocity

        "water", #external fluid
        cold_fluid_inlet_temp, #external temp
        1.0u"bar", #external pressure
        1.0u"m/s", #external fluid velocity

        pipe_diameter * 1.0u"cm", #pipe ID
        (pipe_diameter + 1) * 1.0u"cm", #pipe OD
        1.0u"m", #pipe length
        0.0015u"mm", #pipe absolute roughness
        [copper_layer]
    ) / pipe_length) - Q)

    push!(Q_per_m_vals, Q)
 
    if i == 1
        prev_answer = 1
    else
        prev_answer = ustrip(corresponding_pipe_diameter_vals[i-1])
    end

    push!(corresponding_pipe_diameter_vals, find_zero(find_pipe_diameter_from_heat_duty, prev_answer) * 1.0u"cm")
end

Q_per_m_vals

corresponding_pipe_diameter_vals

using Plots
plot(ustrip.(corresponding_pipe_diameter_vals), Q_per_m_span)
=#

# -- Problem Parameters --
s_values = range(0.0u"m", 1.0u"m", 50)

max_temp_over_start = 50.0u"K"

temp_step = max_temp_over_start / s_values[end]

inner_diameters = []
outer_diameters = []
stored_previous_diameter = 0.0

for (i, s) in enumerate(s_values) 
    if i == 1
        inner_diameters = []
        outer_diameters = []
        stored_previous_diameter = 1.0u"cm"
    end
    
    # -- Wall thickness Calculations 
    internal_pressure = 1.0u"bar"
    joint_efficiency_factor = 0.75
    material_yield_strength = 205u"MPa"

    wall_thickness = (internal_pressure * stored_previous_diameter) / (2.0 * joint_efficiency_factor * (0.8 * material_yield_strength))

    if wall_thickness <= 1.0u"cm"
        wall_thickness = 1.0u"cm"
    end

    stripped_wall_thickness = ustrip(uconvert(u"cm", wall_thickness))

    # -- inner diameter calculations -- 
    temp_at_s = 310.0u"K" + temp_step * s

    desired_Q_per_m = 1000.0u"W/m"

    find_pipe_diameter_from_heat_duty = pipe_diameter -> ((total_heat_duty(
        temp_at_s, #hot side temp
        1.0u"kg/s", #hot side mass flow
        300u"K", #cold side temp
        1.0u"kg/s", #cold side mass flow
        "water", #internal fluid
        temp_at_s, #internal tempm
        1.0u"bar", #internal pressure
        1.0u"m/s", #internal fluid velocity

        "water", #external fluid
        300u"K", #external temp
        1.0u"bar", #external pressure
        1.0u"m/s", #external fluid velocity

        pipe_diameter * 1.0u"cm", #pipe ID
        (pipe_diameter + stripped_wall_thickness) * 1.0u"cm", #pipe OD
        1.0u"m", #pipe length
        0.0015u"mm", #pipe absolute roughness
        [copper_layer]
    ) / 1.0u"m") - desired_Q_per_m)
    
    diameter_required = find_zero(find_pipe_diameter_from_heat_duty, 1.0) #* 1.0u"cm"

    stored_previous_diameter = diameter_required * 1.0u"cm"

    push!(inner_diameters, diameter_required)
    push!(outer_diameters, diameter_required + stripped_wall_thickness)
end

using Interpolations

inner_diameters

outer_diameters

get_inner_diameter_at_s = linear_interpolation(ustrip.(s_values), inner_diameters, extrapolation_bc=Line())

get_outer_diameter_at_s = linear_interpolation(ustrip.(s_values), outer_diameters, extrapolation_bc=Line())

#can accept a value from 0.0 to 1.0 meters
get_inner_diameter_at_s(100)

using LinearAlgebra

# The Procedural Generator
function generate_reactor_sdf(x::Float64, y::Float64, z::Float64)
    #x_clamped = clamp(x, 0.0, 1.0)

    #if x < 0.0 || x > 1.0
        #return 0.1 # Return positive (air) if we are off the ends
    #end

    target_inner_radius = get_inner_diameter_at_s(x) / 2
    target_outer_radius = get_outer_diameter_at_s(x) / 2

    dist_from_axis = sqrt((y - 0.0)^2 + (z - 0.0)^2)
    
    d_inner_cylinder = dist_from_axis - target_inner_radius
    d_outer_cylinder = dist_from_axis - target_outer_radius

    #println("inner cylinder dist, ", d_inner_cylinder, " point dist from axis, ", dist_from_axis, " outer cylinder dist", d_outer_cylinder)

    #=
    if target_inner_radius <= dist_from_axis < target_outer_radius
        return 0.0
    else
        return 1.0
        #return max(d_outer_cylinder, -d_inner_cylinder)
    end
    =#
    
    #just using max(d_outer_cylinder, -d_inner_cylinder) doesn't return anything meaningful
    #it seems like d_inner_cylinder is returning bogus, and min(d_inner_cylinder, d_outer_cylinder) looks the same as d_outer_cylinder
    return min(d_inner_cylinder, d_outer_cylinder)
end

using GLMakie
#=
fig = Figure(size = (800, 600))
ax = Axis3(fig[1, 1], 
    aspect = :data, 
    title = "Generated Reactor Geometry",
    xlabel = "Length (m)", ylabel = "Width (m)", zlabel = "Height (m)"
)

# Use `contour` to draw the isosurface where SDF == 0 (the skin of the pipe)
contour!(ax, x_range[1] .. x_range[end], y_range[1] .. y_range[end], z_range[1] .. z_range[end], vol_data, 
    levels = [0.0],       # The surface is at distance 0
    color = :orange, 
    transparency = false,
    
)
=#

# Optional: Add wireframe to see the curvature changes


x_range = range(0, 100, 100)

x_range = range(-0.1, 0.9, length=100) # Length
y_range = range(-1, 1, length=50)  # Width
z_range = range(-1, 1, length=50)  # Height

vol_data = [generate_reactor_sdf(x, y, z) for x in x_range, y in y_range, z in z_range]

function find_closest_to_zero(matrix)
    # Calculate the absolute values of all elements in the matrix
    abs_matrix = abs.(matrix)

    # Use findmin to get the minimum absolute value and its linear index
    min_abs_value, linear_index = findmin(abs_matrix)

    # Return the original value using the found linear index
    return matrix[linear_index], linear_index
end

count(!iszero, vol_data)

find_closest_to_zero(vol_data)

any(x -> isapprox(x, 0.0), vol_data)

contour(
    x_range, y_range, z_range, generate_reactor_sdf, 
    transparency = true,
    color = (:black, 0.01),
    levels=[0]
)

#display(fig)