using Unitful
using CoolProp
using Roots

function coolprop_mixture(species, mole_fractions)
    total_fraction = sum(mole_fractions)
    corrected_mole_fractions = [mole_fraction / total_fraction for mole_fraction in mole_fractions]

    coolprop_mixture_string = ""
    for i in eachindex(species)
        if i != length(species)
            coolprop_mixture_string *= species[i] * "[" * string(corrected_mole_fractions[i]) * "]" * "&"
        else 
            coolprop_mixture_string *= species[i] * "[" * string(corrected_mole_fractions[i]) * "]"
        end
    end
    return coolprop_mixture_string
end

water_methane_mixture = coolprop_mixture(["water", "methane"], [3.0, 1.0])

function reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, characteristic_length)
    fluid_density = PropsSI("D", "T", fluid_temp, "P", fluid_pressure, fluid)
    dynamic_viscosity = PropsSI("V", "T", fluid_temp, "P", fluid_pressure, fluid)
    return (fluid_density * fluid_velocity * characteristic_length) / (dynamic_viscosity)
end

Re_range = 2200:1:100000

f_vals = []

for Re in Re_range
    find_f = x -> 2 * log10(Re*sqrt(x)) - 0.8 - 1/sqrt(x)
    f = find_zero(find_f, 1e-3)
    push!(f_vals, f)
end
#wow, this only takes 0.66 seconds!

using Interpolations
Re2f = linear_interpolation(Re_range, f_vals)

channel_diameter = 3.0u"mm"

channel_radius = channel_diameter / 2

area = (pi * channel_radius^2) / 2

wetted_perimeter = 1.22 * channel_radius

hydraulic_diameter = (4 * area) / wetted_perimeter

fluid_temp = 700u"°C"

fluid_pressure = 1.0u"bar"

pipe_length = 1.0u"m" 

fluid_density = fluid_density = PropsSI("Dmass", "T", fluid_temp, "P", fluid_pressure, water_methane_mixture)

PropsSI("V", "T", fluid_temp, "P", fluid_pressure, water_methane_mixture)

fluid_velocity = 0.1u"m/s"

Re = reynolds_number(water_methane_mixture, fluid_temp, fluid_pressure, fluid_velocity, hydraulic_diameter)

Re = ustrip(Re)

f = Re2f(Re)

ΔP = f * (pipe_length / hydraulic_diameter) * (fluid_density * fluid_velocity^2 / 2)

ΔP = uconvert(u"Pa", ΔP)

ΔP_per_m = f * (fluid_density * fluid_velocity^2) / (2 * hydraulic_diameter)

ΔP_per_m = uconvert(u"Pa/m", ΔP_per_m)
