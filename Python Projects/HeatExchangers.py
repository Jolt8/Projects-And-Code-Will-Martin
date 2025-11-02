import CoolProp.CoolProp as cp
from CoolProp.CoolProp import PropsSI
import scipy.odr
import scipy.optimize


#from PressureDropCalculations import reynolds_number
from PressureDropCalculations import reynolds_number, cantera_reynolds_number

import ht as ht 
from ht import insulation, conduction, conv_external, conv_internal, boiling_plate

#from CoolProp.CoolProp import IProps, get_Fluid_in dex
#from CoolProp.CoolProp import PhaseSI
from CoolProp.HumidAirProp import HAPropsSI

import colebrook

import numpy as np

import math 

import scipy
from scipy import optimize

from pyfluids import Fluid, FluidsList, Input 

from plots import plot_3d, plot_2d

import cantera as ct

stream = ct.Solution("gri30.yaml", "gri30")
stream.TPX = 300, 100000, "H2O:1"

stream = ct.Water()
stream.TPX = 300, 100000, "H2O:1"

def prandtl_number(fluid, fluid_temperature, fluid_pressure):
    dynamic_viscosity = cp.PropsSI("V", "T", fluid_temperature, "P", fluid_pressure, fluid)
    thermal_cond = cp.PropsSI("conductivity", "T", fluid_temperature, "P", fluid_pressure, fluid)
    specific_heat = cp.PropsSI("C", "T", fluid_temperature, "P", fluid_pressure, fluid)
    density = cp.PropsSI("D", "T", fluid_temperature, "P", fluid_pressure, fluid)
    return (dynamic_viscosity * specific_heat) / (thermal_cond)
#print("prantl", prandtl_number("water", 300, 100000))

def cantera_prandtl_number(solution):
    return (solution.viscosity * solution.cp_mass) / (solution.thermal_conductivity)

def cantera_grashof_number(solution, outer_characteristic_length, surface_temperature):
    """
        Note: you are inputting a cantera solution object that represents ambient conditions, most likely air at 300K and ct.one_atm
    """
    gravitational_acceleration = 9.81
    
    thermal_expansion_coeff = solution.thermal_expansion_coeff
    kinematic_viscosity = solution.viscosity / solution.density
    
    #the abs is to make sure it works for cases where wall temperature is less than ambient fluid temperature
    return (gravitational_acceleration * thermal_expansion_coeff * (abs(surface_temperature - solution.T)) * (outer_characteristic_length ** 3)) / (kinematic_viscosity ** 2)
    
def grashof_number(ambient_fluid, ambient_fluid_temperature, ambient_fluid_pressure, characteristic_length, heat_exchange_surface_temperature):
    gravitational_acceleration = 9.81
    coefficient_of_expansion = PropsSI("isobaric_expansion_coefficient", "T", ambient_fluid_temperature, "P", ambient_fluid_pressure, ambient_fluid)
    
    dynamic_viscosity = cp.PropsSI("V", "T", ambient_fluid_temperature, "P", ambient_fluid_pressure, ambient_fluid)
    density = cp.PropsSI("D", "T", ambient_fluid_temperature, "P", ambient_fluid_pressure, ambient_fluid)
    kinematic_viscosity = dynamic_viscosity / density
    
    #the abs is to make sure it works for cases where wall temperature are lesser than ambient fluid temperature
    return (gravitational_acceleration * coefficient_of_expansion * (abs(heat_exchange_surface_temperature - ambient_fluid_temperature)) * (characteristic_length ** 3)) / (kinematic_viscosity ** 2)
    
#print ("grashof", grashof_number("Water", 320, 100000, 0.1, 350))



def ht_internal_conv(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_characteristic_length, pipe_length, absolute_pipe_roughness):
    """Finds the nusselt number for internal convection (convection from liquid to the pipe)

    Args:
        fluid (str): coolprop compatible fluid
        T (K): fluid temperature
        P (p): fluid pressure
        fluid_velocity (m/s): _description_
        pipe_diameter (m): _description_
        pipe_length (m): _description_
        pipe_roughness (mm): Absolute roughness of the pipe (Copper pipe is around 0.0015mm)

    Returns:
        Nusselt Number (dimensionless)
    """
    fluid_density = cp.PropsSI("D", "T", fluid_temp, "P", fluid_pressure, fluid)
    Re = reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_characteristic_length)
    prandtl = prandtl_number(fluid, fluid_temp, fluid_pressure)
    relative_pipe_roughness = (absolute_pipe_roughness / pipe_characteristic_length)
    friction_factor = colebrook.bntFriction(Re, relative_pipe_roughness, 4)
    return ht.conv_internal.Nu_conv_internal(Re, prandtl, relative_pipe_roughness, pipe_characteristic_length, pipe_length, friction_factor)
#print ("pipe_heat_loss_internal", ht_internal_conv("Water", 300, 100000, 3, 0.0125, 3, 0.0013))


def cantera_ht_internal_conv(solution, fluid_velocity, inner_characteristic_length, pipe_length, absolute_pipe_roughness):
    """Finds the nusselt number for internal convection (convection from liquid to the pipe)

    Args:
        solution: cantera object
        fluid_velocity 
        inner_characteristic_length (float): the inner diameter of the pipe
        pipe_length
        absolute_pipe_roughness (float): reminder, this is not relative pipe roughness, only absolute

    Returns:
        _type_: _description_
    """
    
    
    Re = cantera_reynolds_number(solution, fluid_velocity, inner_characteristic_length)
    prandtl = cantera_prandtl_number(solution)
    relative_pipe_roughness = (absolute_pipe_roughness / inner_characteristic_length)
    friction_factor = colebrook.bntFriction(Re, relative_pipe_roughness, 4)
    return ht.conv_internal.Nu_conv_internal(Re, prandtl, relative_pipe_roughness, inner_characteristic_length, pipe_length, friction_factor)
#print ("pipe_heat_loss_internal", ht_internal_conv("Water", 300, 100000, 3, 0.0125, 3, 0.0013))

def ht_external_pipe(external_fluid, external_fluid_temp, external_fluid_pressure, fluid_velocity, pipe_diameter):
    """Finds the nusselt number for external convection (convection from the pipe to a medium outside of the pipe)
    
    Args:
        fluid (str): coolprop compatible fluid
        T (K): fluid temperature
        P (p): fluid pressure
        fluid_velocity (m/s): _description_
        pipe_diameter (m): _description_
        pipe_length (m): _description_
        pipe_roughness (mm): Absolute roughness of the pipe (Copper pipe is around 0.0015mm)
    
    Returns:
        Nusselt Number (dimensionless)
    """
    Re = reynolds_number(external_fluid, external_fluid_temp, external_fluid_pressure, fluid_velocity, pipe_diameter)
    
    fluid_prandtl = prandtl_number(external_fluid, external_fluid_temp, external_fluid_pressure)
    
    external_prandtl = prandtl_number(external_fluid, external_fluid_temp, external_fluid_pressure)
    fluid_dynamic_viscosity = cp.PropsSI("V", "T", external_fluid_temp, "P", external_fluid_pressure, external_fluid)
    
    external_dynamic_viscosity = cp.PropsSI("V", "T", external_fluid_temp, "P", external_fluid_pressure, external_fluid)
    return ht.conv_external.Nu_external_cylinder(Re, fluid_prandtl, external_prandtl, fluid_dynamic_viscosity, external_dynamic_viscosity)
#print ("pipe_heat_loss_external", ht_external_pipe("Water", 300, 100000, 3, 0.01))


def cantera_ht_external_pipe(solution, inner_characteristic_length, pipe_length, 
    absolute_pipe_roughness, external_fluid_temp_at_wall, fluid_velocity):
    """Finds the nusselt number for external convection (convection from the pipe to a medium outside of the pipe)

    Args:
        solution: cantera solution object
        inner_characteristic_length (m): diameter of pipe
        pipe_length (m):
        absolute_pipe_roughnessfile (mm): Absolute roughness of the pipe (Copper pipe is around 0.0015mm)
        external_fluid_temp (K): temperature of fluid outside of pipe
        external_fluid_temp_at_wall (K):
        external_fluid_pressure (Pa): pressure of fluid outside of pipe
        fluid_velocity (m/s): 
        dynamic_viscosity (_type_): 
        pipe_diameter (_type_): _description_

    Returns:
        nusselt_number: dimensionless
    """
    Re = cantera_reynolds_number(solution, fluid_velocity, inner_characteristic_length)
    
    fluid_prandtl = cantera_prandtl_number(solution)
    
    pre_wall_viscosity = solution.viscosity 
    
    solution.T = external_fluid_temp_at_wall
    wall_viscosity = solution.viscosity
    
    external_prandtl = cantera_prandtl_number(solution)

    return ht.conv_external.Nu_external_cylinder(Re, fluid_prandtl, external_prandtl, pre_wall_viscosity, wall_viscosity)
#print ("pipe_heat_loss_external", ht_external_pipe("Water", 300, 100000, 3, 0.01))

def ht_external_plate(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_characteristic_length):
    """Finds the nusselt number for internal convection (convection from plate heat exchanger to a medium outside of pipe)

    Args:
        fluid (str): coolprop compatible fluid
        T (K): fluid temerpature
        P (p): fluid pressure
        fluid_velocity (m/s): _description_
        pipe_diameter (m): _description_
        pipe_length (m): _description_
        pipe_roughness (mm): Absolute roughness of the pipe (Copper pipe is around 0.0015mm)

    Returns:
        Nusselt Number (dimensionless)
    """
    
    Re = reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_characteristic_length)
    fluid_prandtl = prandtl_number(fluid, fluid_temp, fluid_pressure)
    return ht.conv_external.Nu_external_horizontal_plate(Re, fluid_prandtl, pipe_characteristic_length)
#print("plate_heat_loss_external", ht_external_plate("Water", 300, 100000, 3, 0.01))


def cantera_ht_external_plate(solution, fluid_velocity, pipe_characteristic_length):
    """Finds the nusselt number for internal convection (convection from plate heat exchanger to a medium outside of pipe)

    Args:
        fluid (str): coolprop compatible fluid
        T (K): fluid temerpature
        P (p): fluid pressure
        fluid_velocity (m/s): _description_
        pipe_diameter (m): _description_
        pipe_length (m): _description_
        pipe_roughness (mm): Absolute roughness of the pipe (Copper pipe is around 0.0015mm)

    Returns:
        Nusselt Number (dimensionless)
    """
    
    Re = cantera_reynolds_number(solution, fluid_velocity, pipe_characteristic_length)
    fluid_prandtl = cantera_prandtl_number(solution)
    return ht.conv_external.Nu_external_horizontal_plate(Re, fluid_prandtl, pipe_characteristic_length)
#print("plate_heat_loss_external", ht_external_plate("Water", 300, 100000, 3, 0.01))


def ht_submerged_coil(ambient_fluid, ambient_fluid_temperature, ambient_fluid_pressure, heat_exchange_surface_temperature, outer_characteristic_length, horizontal=False):
    """Finds the nusselt number for pipes submerged in a fluid

    Args:
        fluid (str): coolprop compatible fluid
        T (K): fluid temerpature
        P (p): fluid pressure
        fluid_velocity (m/s): _description_
        pipe_diameter (m): _description_
        pipe_length (m): _description_
        pipe_roughness (mm): Absolute roughness of the pipe (Copper pipe is around 0.0015mm)

    Returns:
        Nusselt Number (dimensionless)
    """
    Gr = grashof_number(ambient_fluid, ambient_fluid_temperature, ambient_fluid_pressure, outer_characteristic_length, heat_exchange_surface_temperature)
    fluid_prandtl = prandtl_number(ambient_fluid, ambient_fluid_temperature, ambient_fluid_pressure)
    return ht.conv_free_immersed.Nu_coil_Xin_Ebadian(fluid_prandtl, Gr, horizontal)
#print("Submerged coil heat loss external", ht_submerged_coil("Water", 350, 100000, 400, 0.125, False))


def cantera_ht_submerged_coil(solution, heat_exchange_surface_temperature, outer_characteristic_length, horizontal=False):
    """Finds the nusselt number for pipes submerged in a fluid

    Args:
        fluid (str): coolprop compatible fluid
        T (K): fluid temerpature
        P (p): fluid pressure
        fluid_velocity (m/s): _description_
        pipe_diameter (m): _description_
        pipe_length (m): _description_
        pipe_roughness (mm): Absolute roughness of the pipe (Copper pipe is around 0.0015mm)

    Returns:
        Nusselt Number (dimensionless)
    """
    Gr = cantera_grashof_number(solution, outer_characteristic_length, heat_exchange_surface_temperature)
    fluid_prandtl = cantera_prandtl_number(solution)
    return ht.conv_free_immersed.Nu_coil_Xin_Ebadian(fluid_prandtl, Gr, horizontal)
#print("Submerged coil heat loss external", ht_submerged_coil("Water", 350, 100000, 400, 0.125, False))

def ht_boiling_flow_chen_bennet(fluid, fluid_phase_temperature, fluid_phase_pressure, mass_flow_rate, 
                                gas_phase_temperature, gas_phase_pressure,
                                pipe_characteristic_length, quality_at_exchanger_specific_interval, 
                                wall_temperature
                                ):
    
    fluid_phase_density = cp.PropsSI("D", "T", fluid_phase_temperature, "P", fluid_phase_pressure, fluid)
    gas_density = cp.PropsSI("D", "T", gas_phase_temperature, "P|gas", gas_phase_pressure, fluid)
    
    fluid_phase_viscosity = cp.PropsSI("V", "T", fluid_phase_temperature, "P", fluid_phase_pressure, fluid)
    gas_phase_viscosity = cp.PropsSI("V", "T", gas_phase_temperature, "P", gas_phase_pressure, fluid)
    
    fluid_phase_thermal_conductivity = cp.PropsSI("conductivity", "T", fluid_phase_temperature, "P", fluid_phase_pressure, fluid)
    fluid_phase_specific_heat_capacity = cp.PropsSI("C", "T", fluid_phase_temperature, "P", fluid_phase_pressure, fluid)
    fluid_phase_heat_of_vaporization = (cp.PropsSI("H", "P", 100000, "Q", 1, fluid) - cp.PropsSI("H", "P", 100000, "Q", 0, fluid)) / 1000
    fluid_surface_tension = cp.PropsSI("surface_tension", "T", fluid_phase_temperature, "Q", 0, fluid)
    
    delta_saturation_pressure = cp.PropsSI("P", "T", gas_phase_temperature, "Q", 1, fluid) - cp.PropsSI("P", "T", fluid_phase_temperature, "Q", 0, fluid)
    
    """
    volumetric_flow_rate = fluid_velocity * ((pipe_characteristic_length ** 2) * math.pi)
    mass_flow_rate = volumetric_flow_rate * fluid_phase_temperature
    """
    
    return ht.boiling_flow.Chen_Bennett(mass_flow_rate, quality_at_exchanger_specific_interval, pipe_characteristic_length, 
                                        fluid_phase_density, gas_density,
                                        fluid_phase_viscosity, gas_phase_viscosity,
                                        fluid_phase_thermal_conductivity, fluid_phase_specific_heat_capacity, fluid_phase_heat_of_vaporization,
                                        fluid_surface_tension, delta_saturation_pressure, 
                                        wall_temperature
                                        )
#print("boiling flow", ht_boiling_flow_chen_bennet(
    fluid="water",
    fluid_phase_temperature=300,
    fluid_phase_pressure=100000,
    mass_flow_rate=0.5,
    gas_phase_temperature=450,
    gas_phase_pressure=300000,
    pipe_characteristic_length=0.1,
    quality_at_exchanger_specific_interval=0.4,
    wall_temperature=3,
#))

def cantera_ht_boiling_flow_chen_bennet(condensed_solution, boiled_solution, quality_at_exchanger_specific_interval, fluid_conductivity, mass_flow_rate, 
                                inner_characteristic_length, 
                                wall_temperature
                                ):
    """
    Notes: 
    - requires two solution objects representing the condensed and boiled state 
    - the solutions used must have EOS and transport properties (ex. thermal_conductivity, quality, etc.)
    - the solution = ct.water() is an example of a fluid that has all the necessary parameters
    
    
    """
    #psat at wall temperature - psat at internal fluid temperature
    
    """
    volumetric_flow_rate = fluid_velocity * ((pipe_characteristic_length ** 2) * math.pi)
    mass_flow_rate = volumetric_flow_rate * fluid_phase_temperature
    """
    
    return ht.boiling_flow.Chen_Bennett(mass_flow_rate, quality_at_exchanger_specific_interval, inner_characteristic_length, 
                                        condensed_solution.density, boiled_solution.density,
                                        condensed_solution.viscosity, boiled_solution.density,
                                        condensed_solution.thermal_conductivity, condensed_solution.cp_mass, condensed_solution.hvap,
                                        #FIXME: I'm not sure that cantera has a surface tension attribute
                                        condensed_solution.surface_tension, abs(condensed_solution.P_sat - boiled_solution.P_sat), 
                                        wall_temperature
                                        )
#print("boiling flow", ht_boiling_flow_chen_bennet(
    fluid="water",
    fluid_phase_temperature=300,
    fluid_phase_pressure=100000,
    mass_flow_rate=0.5,
    gas_phase_temperature=450,
    gas_phase_pressure=300000,
    pipe_characteristic_length=0.1,
    quality_at_exchanger_specific_interval=0.4,
    wall_temperature=3,
#))



def heat_transfer_coeff(Nu, fluid, fluid_temp, fluid_pressure, characteristic_length):
    fluid_thermal_cond = cp.PropsSI("conductivity", "T", fluid_temp, "P", fluid_pressure, fluid)
    return (Nu * fluid_thermal_cond) / characteristic_length
#print("single heat transfer coefficient:", heat_transfer_coeff(ht_external_pipe("Water", 300, 100000, 3, 0.01), "Water", 300, 100000, 0.01))


def cantera_heat_transfer_coeff(solution, inner_characteristic_length, Nu):
    return (Nu * solution.thermal_conductivity) / inner_characteristic_length

def overall_heat_transfer_coeff(internal_nusselt_number, internal_fluid, internal_fluid_temp, internal_fluid_pressure, internal_characteristic_length, internal_area, 
                                        external_nusselt_number, external_fluid, external_fluid_temp, external_fluid_pressure, external_characteristic_length, external_area, 
                                        layers):
    """
    Keep in mind that layers if formatted [(thickness, thermal conductivity, area), (thickness, thermal conductivity, area)]
    """
    conductive_resistance = sum(thickness / (conductivity * area) for thickness, conductivity, area in layers)
    
    internal_heat_transfer_coeff = heat_transfer_coeff(internal_nusselt_number, internal_fluid, internal_fluid_temp, internal_fluid_pressure, internal_characteristic_length)
    external_heat_transfer_coeff = heat_transfer_coeff(external_nusselt_number, external_fluid, external_fluid_temp, external_fluid_pressure, external_characteristic_length)
    
    return 1 / ((1 / (internal_heat_transfer_coeff * internal_area)) + (conductive_resistance) + (1 / (external_heat_transfer_coeff * external_area)))

def cantera_overall_heat_transfer_coeff(internal_nusselt_number, internal_solution, inner_characteristic_length, internal_area, 
                                        external_nusselt_number, external_solution, outer_characteristic_length, external_area, 
                                        layers):
    """
    Keep in mind that layers if formatted [(thickness 1, thermal conductivity 1, area 1), (thickness 2, thermal conductivity 2, area 2)]
    the thickness, thermal conductivity, and area are all for the materials seperating the two fluids
    """
    conductive_resistance = sum(thickness / (conductivity * area) for thickness, conductivity, area in layers)
    
    internal_heat_transfer_coeff = cantera_heat_transfer_coeff(internal_solution, inner_characteristic_length, internal_nusselt_number)
    external_heat_transfer_coeff = cantera_heat_transfer_coeff(external_solution, outer_characteristic_length, external_nusselt_number)
    
    return 1 / ((1 / (internal_heat_transfer_coeff * internal_area)) + (conductive_resistance) + (1 / (external_heat_transfer_coeff * external_area)))

#print("Overall heat transfer v2", overall_heat_transfer_coeff(ht_internal_conv("Isobutane", 300, 500000, 23, 0.125, 3, 0.0015), "Isobutane", 300, 500000, 0.125, 0.1, 
                                                                  #ht_external_pipe("Water", 300, 100000, 3, 0.3), "water", 300, 100000, 0.3, 0.3, 
                                                                  #[(0.012, 400, 0.125)]))
#print("Overall heat transfer v2", overall_heat_transfer_coeff(
    internal_nusselt_number= ht_internal_conv("Isobutane", 300, 500000, 23, 0.0125, 3, 0.0015),
    internal_fluid="Isobutane", 
    internal_fluid_temp=(10 + 273), 
    internal_fluid_pressure=500000, 
    internal_characteristic_length=0.125, 
    internal_area=0.1,
    external_nusselt_number= ht_submerged_coil("Water", 400, 100000, 350, 0.009525, False),
    external_fluid="Water", 
    external_fluid_temp=(90 + 273), 
    external_fluid_pressure=100000, 
    external_characteristic_length=0.009525, 
    external_area=0.03, 
    layers=[(0.0016, 400, 0.01)], 
#))
"""
print("Overall heat transfer v3", overall_heat_transfer_coeff(
    internal_nusselt_number= ht_internal_conv("Isobutane", 300, 500000, 23, 0.0125, 3, 0.0015),
    internal_fluid="Isobutane", 
    internal_fluid_temp=300, 
    internal_fluid_pressure=500000, 
    internal_characteristic_length=0.125, 
    internal_area=0.1,
    external_nusselt_number= ht_submerged_coil("Water", 350, 100000, 300, 0.0125, False),
    external_fluid="Water", 
    external_fluid_temp=300, 
    external_fluid_pressure=100000, 
    external_characteristic_length=0.0125, 
    external_area=0.1, 
    layers=[(0.0016, 400, 0.01)]
#))

#print("Overall heat transfer v4", overall_heat_transfer_coeff(
    internal_nusselt_number= ht_boiling_flow_chen_bennet(
        fluid="Water",
        fluid_phase_temperature=300,
        fluid_phase_pressure=100000,
        mass_flow_rate=0.5,
        gas_phase_temperature=390,
        gas_phase_pressure=300000,
        pipe_characteristic_length=0.0125,
        quality_at_exchanger_specific_interval=0.4,
        wall_temperature=3,
        ),
    internal_fluid="Water", 
    internal_fluid_temp=300, 
    internal_fluid_pressure=500000, 
    internal_characteristic_length=0.0125, 
    internal_area=0.1,
    external_nusselt_number= ht_submerged_coil(
        ambient_fluid="Water",
        ambient_fluid_temperature=400,
        ambient_fluid_pressure=100000,
        heat_exchange_surface_temperature=300,
        outer_characteristic_length=0.0125,
        horizontal=False,
    ),
    external_fluid="Water", 
    external_fluid_temp=300, 
    external_fluid_pressure=100000, 
    external_characteristic_length=0.0125, 
    external_area=0.1, 
    layers=[(0.0016, 400, 0.01)]
#))
                                                         
#300, 100000, 3, 1, 0.01, 0.0013, 400, 0.01))
"""

""" Functions needed for Pipe heat loss calculations
 - Internal convection (convection from liquid to pipe wall)
 - External Convection (convection on the outside of the pipe)
 - Conduction and Shape Factors
 - Heat Transfer By Radiation
 - Database for Insulating and Refractory metals



 Functions needed for Heat Exchanger Calculations
 - Heat Exchanger Sizing and Rating
 - Convection to Plate Heat Exchangers 
 - Conduction and Shape Factors
 - Free convection to Enclosed Bodies
 - Heat Transfer and Pressure Drop across tube bundles
 - Condensation (Only if condensation takes place on the walls)
 - Boiling in plate and frame exchangers (ht.boiling_plate)
"""

def log_mean_temperature_difference(T_hot_in, T_hot_out, T_cold_in, T_cold_out):
    T1 = abs(T_hot_in - T_cold_in)
    T2 = abs(T_hot_out - T_cold_out)
    return (T1 - T2) / (np.log(T1 / T2))
#print ("ln_temp", log_mean_temperature_difference(90 + 273, 70 + 273, 20 + 273, 60 + 273))



def heat_duty(heat_transfer_coefficient, inner_characteristic_length, pipe_length, T_hot_in, T_hot_out, T_cold_in, T_cold_out):
    area = (math.pi * inner_characteristic_length) * pipe_length
    log_mean_temp_diff = log_mean_temperature_difference(T_hot_in, T_hot_out, T_cold_in, T_cold_out)
    return heat_transfer_coefficient * log_mean_temp_diff

def cantera_air_heat_loss(internal_solution, expected_temperature_change, 
                  mass_flow,
                  outer_wall_temperature,
                  inner_characteristic_length, outer_characteristic_length, pipe_length, absolute_pipe_roughness, pipe_conductivity,
                  air_temperature_initial, air_temperature_final):
    """Gets the wattage of heat transfer for the parameters given (usually used for environmental heat loss)

    Args:
        internal_solution (obj): cantera solution object
        internal_temperature_initial (K): initial temperature of the fluid object, you can usually just make this solution.T
        internal_temperature_final (K): final temperature of the fluid object, make it solution.T + expected temperature rise
        mass_flow (kg/s)
        inner_wall_temperature (K): expected temperature of the inside wall
        outer_wall_temperature (K): expected temeperature of outside wall
        inner_characteristic_length (m): pipe inside diameter
        outer_characteristic_length (m): pipe outside diameter
        pipe_length (m)
        absolute_pipe_roughness (mm): remember, this is not the same as relative pipe roughness
        pipe_conductivity (W/m*K): conductivity of the pipe
        air_temperature_initial (_type_): initial tempertaure of air
        air_temperature_final (_type_): expected final temperature of air

    Returns:
        heat transfer to environment in watts
        Note: to find the enthalpy change, multiply the watts by the time spent in the pipe (pipe length (m) / fluid velocity (m/s))
    """
    volumetric_flow = mass_flow / internal_solution.density
    area = math.pi * ((inner_characteristic_length / 2) ** 2)
    velocity = volumetric_flow / area
    air_solution = ct.Solution("air.yaml")
    air_solution.TPX = air_temperature_initial, ct.one_atm, "N2:0.78, O2:0.21"
    print(cantera_ht_internal_conv(internal_solution, velocity, inner_characteristic_length, pipe_length, absolute_pipe_roughness))
    print(cantera_ht_submerged_coil(air_solution, outer_wall_temperature, outer_characteristic_length, False))
    print(internal_solution.T, internal_solution.T + expected_temperature_change, air_temperature_initial, air_temperature_final)
    print(cantera_overall_heat_transfer_coeff(
        internal_nusselt_number=cantera_ht_internal_conv(internal_solution, velocity, inner_characteristic_length, pipe_length, absolute_pipe_roughness),
        internal_solution=internal_solution,
        inner_characteristic_length=inner_characteristic_length, 
        internal_area=(math.pi * inner_characteristic_length * pipe_length),
        #FIXME: I'm not sure what the best way is to get the surface temperature of the heat exchanger 
            #maybe there's a way to estimate it using thermal conductivity, thickness, and LMTD
        external_nusselt_number=cantera_ht_submerged_coil(air_solution, outer_wall_temperature, outer_characteristic_length, False),
        external_solution=air_solution,
        outer_characteristic_length=outer_characteristic_length, 
        external_area=(math.pi * outer_characteristic_length * pipe_length), 
        layers=[((outer_characteristic_length - inner_characteristic_length), pipe_conductivity, (math.pi * outer_characteristic_length * pipe_length))]))
    return heat_duty(
    heat_transfer_coefficient=cantera_overall_heat_transfer_coeff(
        internal_nusselt_number=cantera_ht_internal_conv(internal_solution, velocity, inner_characteristic_length, pipe_length, absolute_pipe_roughness),
        internal_solution=internal_solution,
        inner_characteristic_length=inner_characteristic_length, 
        internal_area=((math.pi * inner_characteristic_length) * pipe_length),
        #FIXME: I'm not sure what the best way is to get the surface temperature of the heat exchanger 
            #maybe there's a way to estimate it using thermal conductivity, thickness, and LMTD
        external_nusselt_number=cantera_ht_submerged_coil(air_solution, outer_wall_temperature, outer_characteristic_length, False),
        external_solution=air_solution,
        outer_characteristic_length=outer_characteristic_length, 
        external_area=((math.pi * outer_characteristic_length) * pipe_length), 
        layers=[((outer_characteristic_length - inner_characteristic_length), pipe_conductivity, ((math.pi * outer_characteristic_length) * pipe_length))]), 
    inner_characteristic_length=inner_characteristic_length,
    pipe_length=pipe_length,
    #area=((math.pi * outer_characteristic_length) * pipe_length),
    T_hot_in=internal_solution.T,
    T_hot_out=internal_solution.T + expected_temperature_change,
    T_cold_in=air_temperature_initial,
    T_cold_out=air_temperature_final
    )

water = ct.Water()
# Set to ~60°C, 1 atm
water.TP = 333.15, ct.one_atm
# internal_solution = water
print(cantera_air_heat_loss(
    internal_solution=water,                     # Cantera Solution for liquid water
    expected_temperature_change=20,               # K — Air temperature rises by 2 K
    mass_flow=0.1,                               # kg/s — Reasonable for air in a duct
    outer_wall_temperature=350,                  # K — Hot water pipe (~77 °C)
    inner_characteristic_length=0.02,            # m — 2 cm pipe inner diameter
    outer_characteristic_length=0.021,            # m — Also 2 cm outer diameter
    pipe_length=1.0,                             # m — 1 meter pipe
    absolute_pipe_roughness=0.000045,            # m — Commercial steel pipe
    pipe_conductivity=50,                        # W/m·K — Stainless steel
    air_temperature_initial=300,                 # K — Ambient air (27 °C)
    air_temperature_final=310                    # K — Small rise due to heat gain
))

print("Boiler", heat_duty(
    overall_heat_transfer_coeff(
        internal_nusselt_number= ht_internal_conv("Water", 30 + 273, 500000, 23, 0.125, 3, 0.0015),
        internal_fluid="Isobutane", 
        internal_fluid_temp=10 + 273, 
        internal_fluid_pressure=500000, 
        internal_characteristic_length=0.125, 
        internal_area=0.1,
        external_nusselt_number= ht_submerged_coil("Water", 100 + 273, 100000, 350, 0.009525, False),
        external_fluid="Water", 
        external_fluid_temp=90 + 273, 
        external_fluid_pressure=100000, 
        external_characteristic_length=0.009525, 
        external_area=0.1, 
        layers=[(0.012, 400, 0.009525)]), 
    inner_characteristic_length=0.01,
    pipe_length=0.5,
    T_hot_in=100 + 273,
    T_hot_out=80 + 273,
    T_cold_in=20 + 273,
    T_cold_out=70 + 273
))
print("Condenser", heat_duty(
    overall_heat_transfer_coeff(
        internal_nusselt_number= ht_internal_conv("Water", 80 + 273, 500000, 23, 0.125, 3, 0.0015),
        internal_fluid="Isobutane", 
        internal_fluid_temp=300, 
        internal_fluid_pressure=500000, 
        internal_characteristic_length=0.125, 
        internal_area=0.1,
        external_nusselt_number= ht_submerged_coil("Water", 10 + 273, 100000, 30 + 273, 0.125, False),
        external_fluid="Water", 
        external_fluid_temp=10 + 273, 
        external_fluid_pressure=100000, 
        external_characteristic_length=0.009, 
        external_area=0.1, 
        layers=[(0.012, 400, 0.125)]), 
    inner_characteristic_length=0.01,
    pipe_length=0.5,
    T_hot_in=90 + 273,
    T_hot_out=50 + 273,
    T_cold_in=10 + 273,
    T_cold_out=30 + 273
))
