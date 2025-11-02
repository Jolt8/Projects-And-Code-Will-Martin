using Clapeyron
using Unitful

model = PR(["Methanol", "Water", "Carbon monoxide", "Carbon Dioxide", "Hydrogen"])

mass_density(model, 1.0u"bar", 573.28u"K", [1.0, 1.0, 1.0, 1.0, 1.0])

bed_void_fraction = 0.3
gas_viscosity = 2.0e-5u"Pa*s"
superficial_mass_velocity = 0.12409u"kg/(m^2*s)"
gas_density = 0.521u"kg/m^3"
catalyst_particle_diameter = 5u"mm"
bed_density = 1000u"kg/m^3"
reactor_cross_sectional_area = 1u"m^2"

term1 = 150 * (1 - bed_void_fraction)^2 / bed_void_fraction^3 * gas_viscosity * superficial_mass_velocity / (gas_density * catalyst_particle_diameter^2)
term1 = uconvert(u"bar/m", term1)
term2 = (1.75 * (1 - bed_void_fraction) * superficial_mass_velocity^2) / (catalyst_particle_diameter * gas_density * bed_void_fraction^3)
term2 = uconvert(u"bar/m", term2)
println("num, ", (1.75 * (1 - bed_void_fraction) * superficial_mass_velocity^2))
println("denom, ", uconvert(u"kg/m^2", (catalyst_particle_diameter * gas_density * bed_void_fraction^3)))
dP_dL = -(term1 + term2) # This is in Pa/m
uconvert(u"bar/m", dP_dL)

bed_density = (catalyst_density * (1 - bed_void_fraction))

dP_dW = dP_dL / (bed_density * reactor_cross_sectional_area) # Pa/kg
uconvert(u"bar/kg", dP_dW)