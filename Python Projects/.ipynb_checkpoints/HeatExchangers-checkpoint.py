import CoolProp.CoolProp as cp
from CoolProp.CoolProp import PropsSI
import scipy.odr
import scipy.optimize


#from PressureDropCalculations import reynolds_number
from PressureDropCalculations import reynolds_number

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

def prandtl_number (fluid, fluid_temperature, fluid_pressure):
    dynamic_viscosity = cp.PropsSI("V", "T", fluid_temperature, "P", fluid_pressure, fluid)
    thermal_cond = cp.PropsSI("conductivity", "T", fluid_temperature, "P", fluid_pressure, fluid)
    specific_heat = cp.PropsSI("C", "T", fluid_temperature, "P", fluid_pressure, fluid)
    density = cp.PropsSI("D", "T", fluid_temperature, "P", fluid_pressure, fluid)
    return (dynamic_viscosity * specific_heat) / (thermal_cond)
print("prantl", prandtl_number("water", 300, 100000))

def grashof_number (ambient_fluid, ambient_fluid_temperature, ambient_fluid_pressure, characteristic_length, heat_exchange_surface_temperature):
    gravitational_acceleration = 9.81
    coefficient_of_expansion = PropsSI("isobaric_expansion_coefficient", "T", ambient_fluid_temperature, "P", ambient_fluid_pressure, ambient_fluid)
    
    dynamic_viscosity = cp.PropsSI("V", "T", ambient_fluid_temperature, "P", ambient_fluid_pressure, ambient_fluid)
    density = cp.PropsSI("D", "T", ambient_fluid_temperature, "P", ambient_fluid_pressure, ambient_fluid)
    kinematic_viscosity = dynamic_viscosity / density
    
    #the abs is to make sure it works for cases where wall temperature are lesser than ambient fluid temperature
    return (gravitational_acceleration * coefficient_of_expansion * (abs(heat_exchange_surface_temperature - ambient_fluid_temperature)) * (characteristic_length ** 3)) / (kinematic_viscosity ** 2)
    
print ("grashof", grashof_number("Water", 320, 100000, 0.1, 350))



def ht_internal_conv(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_diameter, pipe_length, absolute_pipe_roughness):
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
    Re = reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_diameter)
    prandtl = prandtl_number(fluid, fluid_temp, fluid_pressure)
    relative_pipe_roughness = (absolute_pipe_roughness / pipe_diameter)
    friction_factor = colebrook.bntFriction(Re, relative_pipe_roughness, 4)
    return ht.conv_internal.Nu_conv_internal(Re, prandtl, relative_pipe_roughness, pipe_diameter, pipe_length, friction_factor)
print ("pipe_heat_loss_internal", ht_internal_conv("Water", 300, 100000, 3, 0.01, 1, 0.0013))

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
print ("pipe_heat_loss_external", ht_external_pipe("Water", 300, 100000, 3, 0.01))

def ht_external_plate(fluid, fluid_temp, fluid_pressure, fluid_velocity, plate_length):
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
    Re = reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, plate_length)
    fluid_prandtl = prandtl_number(fluid, fluid_temp, fluid_pressure)
    return ht.conv_external.Nu_external_horizontal_plate(Re, fluid_prandtl, plate_length)
print("plate_heat_loss_external", ht_external_plate("Water", 300, 100000, 3, 0.01))


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
    print(Gr)
    fluid_prandtl = prandtl_number(ambient_fluid, ambient_fluid_temperature, ambient_fluid_pressure)
    print(fluid_prandtl)
    return ht.conv_free_immersed.Nu_coil_Xin_Ebadian(fluid_prandtl, Gr, horizontal)
print("Submerged coil heat loss external", ht_submerged_coil("Water", 350, 100000, 400, 0.125, False))



def heat_transfer_coeff(Nu, fluid, fluid_temp, fluid_pressure, characteristic_length):
    fluid_thermal_cond = cp.PropsSI("conductivity", "T", fluid_temp, "P", fluid_pressure, fluid)
    return (Nu * fluid_thermal_cond) / characteristic_length
print("single heat transfer coefficient:", heat_transfer_coeff(ht_external_pipe("Water", 300, 100000, 3, 0.01), "Water", 300, 100000, 0.01))

def overall_heat_transfer_coeff(internal_nusselt_number, internal_fluid, internal_fluid_temp, internal_fluid_pressure, internal_characteristic_length, internal_area, 
                                        external_nusselt_number, external_fluid, external_fluid_temp, external_fluid_pressure, external_characteristic_length, external_area, 
                                        layers):
    """
    Keep in mind that layers if formatted [(thickness 1, thermal conductivity 1, area 1), (thickness 2, thermal conductivity 2, area 2)]
    """
    conductive_resistance = sum(thickness / (conductivity * area) for thickness, conductivity, area in layers)
    
    internal_heat_transfer_coeff = heat_transfer_coeff(internal_nusselt_number, internal_fluid, internal_fluid_temp, internal_fluid_pressure, internal_characteristic_length)
    external_heat_transfer_coeff = heat_transfer_coeff(external_nusselt_number, external_fluid, external_fluid_temp, external_fluid_pressure, external_characteristic_length)
    
    return 1 / ((1 / (internal_heat_transfer_coeff * internal_area)) + (conductive_resistance) + (1 / (external_heat_transfer_coeff * external_area)))

print("Overall heat transfer v2", overall_heat_transfer_coeff(ht_internal_conv("Isobutane", 300, 500000, 23, 0.125, 3, 0.0015), "Isobutane", 300, 500000, 0.125, 0.1, 
                                                                  ht_external_pipe("Water", 300, 100000, 3, 0.3), "water", 300, 100000, 0.3, 0.3, 
                                                                  [(0.012, 400, 0.125)]))

print("Overall heat transfer v2", overall_heat_transfer_coeff(
    internal_nusselt_number= ht_internal_conv("Isobutane", 300, 500000, 23, 0.125, 3, 0.0015),
    internal_fluid="Isobutane", 
    internal_fluid_temp=300, 
    internal_fluid_pressure=500000, 
    internal_characteristic_length=0.125, 
    internal_area=0.1,
    external_nusselt_number=ht_external_pipe("Water", 350, 100000, 3, 0.125),
    external_fluid="Water", 
    external_fluid_temp=300, 
    external_fluid_pressure=100000, 
    external_characteristic_length=0.3, 
    external_area=0.3, 
    layers=[(0.012, 400, 0.125)]
))

print("Overall heat transfer v3", overall_heat_transfer_coeff(
    internal_nusselt_number= ht_internal_conv("Isobutane", 300, 500000, 23, 0.125, 3, 0.0015),
    internal_fluid="Isobutane", 
    internal_fluid_temp=300, 
    internal_fluid_pressure=500000, 
    internal_characteristic_length=0.125, 
    internal_area=0.1,
    external_nusselt_number=ht_submerged_coil("Water", 350, 100000, 400, 0.125, False),
    external_fluid="Water", 
    external_fluid_temp=300, 
    external_fluid_pressure=100000, 
    external_characteristic_length=0.3, 
    external_area=0.3, 
    layers=[(0.012, 400, 0.125)]
))
                                                         
#300, 100000, 3, 1, 0.01, 0.0013, 400, 0.01))

""" Functions needed for Pipe heat loss calculations
 - Internal convection (convection from liquid to pipe wall)
 - External Convection (convection on the outside of the pipe)
 - Conduction and Shape Factors
 - Heat Transfer By Radiation
 - Database for Insulating and Refractory metals

"""

""" Functions needed for Heat Exchanger Calculations
 - Heat Exchanger Sizing and Rating
 - Convection to Plate Heat Exchangers 
 - Conduction and Shape Factors
 - Free convection to Enclosed Bodies
 - Heat Transfer and Pressure Drop across tube bundles
 - Condensation (Only if condensation takes place on the walls)
 - Boiling in plate and frame exchangers (ht.boiling_plate)
"""

def log_mean_temperature_difference(T_hot_in, T_hot_out, T_cold_in, T_cold_out):
    T1 = T_hot_in - T_cold_in
    T2 = T_hot_out - T_cold_out
    return (T1 - T2) / (np.log(T1 / T2))
print ("ln_temp", log_mean_temperature_difference(90 + 273, 70 + 273, 20 + 273, 60 + 273))



def heat_duty(heat_transfer_coefficient, area, T_hot_in, T_hot_out, T_cold_in, T_cold_out):
    log_mean_temp_diff = log_mean_temperature_difference(T_hot_in, T_hot_out, T_cold_in, T_cold_out)
    return heat_transfer_coefficient * area * log_mean_temp_diff

print("useful calculation! v2", heat_duty(500, 0.1, 90 + 273, 70 + 273, 20 + 273, 60 + 273))


#Find this fucntion in ht.py:
#ht.hx.P_NTU_method(m1, m2, Cp1, Cp2, UA=None, T1i=None, T1o=None, T2i=None, T2o=None, subtype='crossflow', Ntp=1, optimal=True) 


