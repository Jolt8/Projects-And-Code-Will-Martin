using Unitful

import CoolProp as Cp

desired_electrical_power = 10u"MW"

cycle_thermal_efficiency = 0.42

reactor_heat_output = desired_electrical_power / cycle_thermal_efficiency

reactor_power_density = 80u"kW/L" #around 50-150 kW/L for a heat pipe reactor 

required_core_volume = uconvert(u"L", reactor_heat_output / reactor_power_density)

reactor_height_to_diameter_ratio = 1

reactor_radius = uconvert(u"m", cbrt(required_core_volume / (2 * pi * reactor_height_to_diameter_ratio)))

reactor_diameter = reactor_radius * 2

reactor_height = reactor_diameter * reactor_height_to_diameter_ratio



