using Unitful
using CoolProp
using Symbolics, Nemo, Groebner

function get_heat_gradient(x)
    return (150u"K/m" * x + 300u"K") 
end
#gradient of 450 at start and 300 at end 

get_heat_gradient(1u"m")

function reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, characteristic_length)
    fluid_density = PropsSI("D", "T", fluid_temp, "P", fluid_pressure, fluid)
    dynamic_viscosity = PropsSI("V", "T", fluid_temp, "P", fluid_pressure, fluid)
    return (fluid_density * fluid_velocity * characteristic_length) / (dynamic_viscosity)
end

# --- TOOLING ---
@variables Re, fluid_density, fluid_velocity, pipe_diameter, dynamic_viscosity

eq = Re ~ (fluid_density * fluid_velocity * pipe_diameter) / (dynamic_viscosity)

try 
    symbolic_linear_solve(eq, pipe_diameter)[1]
catch
    symbolic_solve(eq, pipe_diameter)[1]
end
# --- TOOLING ---

function pipe_diameter_from_reynolds_number(Re, fluid_density, fluid_velocity, dynamic_viscosity)
    return (Re*dynamic_viscosity) / (fluid_density*fluid_velocity)
end

    
function prandtl_number(fluid, fluid_temperature, fluid_pressure)
    dynamic_viscosity = PropsSI("V", "T", fluid_temperature, "P", fluid_pressure, fluid)
    thermal_cond = PropsSI("conductivity", "T", fluid_temperature, "P", fluid_pressure, fluid)
    specific_heat = PropsSI("C", "T", fluid_temperature, "P", fluid_pressure, fluid)
    density = PropsSI("D", "T", fluid_temperature, "P", fluid_pressure, fluid)
    return (dynamic_viscosity * specific_heat) / (thermal_cond)
end

prandtl_number("water", 300u"K", 1u"bar")

function laminar_entry_thermal_Hausen(Re, Pr, pipe_length, pipe_diameter)
    Gz = (pipe_diameter / pipe_length) * Re * Pr
    return 3.66 + (0.0668*Gz)/(1+0.04*(Gz)^(2/3.0))
end

# --- TOOLING ---
@variables Nu, Re, Pr, pipe_length, pipe_diameter

Gz = (pipe_diameter / pipe_length) * Re * Pr

eq = Nu ~ 3.66 + (0.0668*Gz)/(1+0.04*(Gz)^(2/3))

eq_simplified = simplify(eq)

try 
    symbolic_linear_solve(eq_simplified, pipe_diameter)[1]
catch
    symbolic_solve(eq_simplified, pipe_diameter)[1]
end
#I'm starting to realize that depending on this to try an generate new functions is impossible because pipe diameter is present multiple times throughout the equation
# --- TOOLING ---

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
    
    #Keep in mind that layers if formatted [(thickness, thermal conductivity, area), (thickness, thermal conductivity, area)]
    overall_heat_trans_coeff = overall_heat_transfer_coeff(
        internal_NU, internal_thermal_cond, pipe_inner_diameter, internal_area, 
        external_NU, external_thermal_cond, pipe_outer_diameter, external_area, 
        layers
    )
    return overall_heat_trans_coeff
end

struct Layer
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

get_overall_heat_transfer_coefficent(
    "water", #internal fluid
    300u"K", #internal temp
    1u"bar", #internal pressure
    1u"m/s", #internal fluid velocity

    "water", #external fluid
    450u"K", #external temp
    1u"bar", #external pressure
    1u"m/s", #external fluid velocity

    10.0u"cm", #pipe ID
    11.0u"cm", #pipe OD
    1u"m", #pipe length
    0.0015u"mm", #pipe absolute roughness
    [copper_layer]
)






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

function effectiveness_from_NTU(NTU, heat_cap_ratio)
    return (1.0 - exp(-NTU*(1.0 - heat_cap_ratio)))/(1.0 - heat_cap_ratio*exp(-NTU*(1.0 - heat_cap_ratio)))
end

function effectiveness_NTU_method(
    hot_inlet_temp, cp_hot, m_hot, 
    cold_inlet_temp, cp_cold, m_cold, 
    UA)

    cp_min = min_heat_cap_rate(cp_hot, m_hot, cp_cold, m_cold)
    
    NTU = UA / cp_min

    heat_cap_ratio = heat_cap_ratio(cp_hot, m_hot, cp_cold, m_cold)

    effectiveness = effectiveness_from_NTU(NTU, heat_cap_ratio)

    Q = effectiveness * cp_min * (hot_inlet_temp - cold_inlet_temp)
end
