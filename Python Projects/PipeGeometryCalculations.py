from CoolProp.CoolProp import PropsSI
#from CoolProp.CoolProp import IProps, get_Fluid_index
#from CoolProp.CoolProp import PhaseSI
import CoolProp.CoolProp as cp
from CoolProp.HumidAirProp import HAPropsSI

import numpy as np

import math 

import scipy
from scipy import optimize

from pyfluids import Fluid, FluidsList, Input 

def pipe_wall_V_d_outer_d_inner_L(*args):
    def pre_pipe_wall_V_d_outer_d_inner_L(wall_volume, D_outer, D_inner, length):
        """Finds the volume of the pipe walls to determine how much energy it takes to heat up the walls

        Args:
            volume (m3): volume of the pipe walls
            D_outer (m): outer diameter of the pipe
            D_inner (m): inner diameter of the pipe
            length (m): total length of the pipe
            
        intermediates:
            out_circle (m): finds the area of the circle for the outer diameter of the pipe
            out_circle (m): finds the area of the circle for the inner diameter of the pipe  

        Returns:
            returns any value with "None"
        """
        out_circle = math.pi*(D_outer/2)**2
        in_circle = math.pi*(D_inner/2)**2
        return ((out_circle - in_circle) * length) - wall_volume
    i = args.index(None)
    return scipy.optimize.fsolve(lambda x: pre_pipe_wall_V_d_outer_d_inner_L(*args[:i], x, *args[i+1:]), 0.1)
print ("pipe wall volume:", pipe_wall_V_d_outer_d_inner_L(None, 0.2, 0.005, 1))

def pipe_inside_V_DL(*args):
    def pre_pipe_inside_V_DL(volume, diameter, length):
        """Finds the volume of the inside of the pipe 

        Args:
            volume (m3): inner volume of the pipe
            diameter (m): inner diameter of the pipe
            length (m): length of the pipe

        Returns:
            returns any value with "None"
        """
        return (math.pi*(diameter/2)**2 * length) - volume
    i = args.index(None)
    return scipy.optimize.fsolve(lambda x: pre_pipe_inside_V_DL(*args[:i], x, *args[i+1:]), 20)

print ("pipe inside volume:", pipe_inside_V_DL(None, 0.01, 7.62))

def pipe_area_A_DL(*args):
    def pre_pipe_area_A_DL(area, diameter, length):
        """Finds the surface area of a pipe for the inner or outer diameter

        Args:
            area (m2): surface area of either the inside or outside of the pipe
            diameter (m): inner or outer diameter of the pipe
            length (m): length of the pipe

        Returns:
            Any value with "None"
        """
        return math.pi*diameter * length - area
    i = args.index(None)
    return scipy.optimize.fsolve(lambda x: pre_pipe_area_A_DL(*args[:i], x, *args[i+1:]), 20)

print ("pipe area:", pipe_area_A_DL(None, 0.0127, 30))


def fluid_velocity_from_litres_second(pipe_inner_diameter, volumetric_flow_rate):
   cross_sectional_area = ((pipe_inner_diameter / 2) ** 2) * math.pi 
   return volumetric_flow_rate * 0.001 / cross_sectional_area

print (fluid_velocity_from_litres_second(0.0127, 3))

def fluid_velocity_from_mass_flow_rate(fluid, fluid_mass_flow_rate, fluid_temperature, fluid_pressure, pipe_inner_diameter):
   density =  cp.PropsSI("D", "T", fluid_temperature, "P", fluid_pressure, fluid)
   volume = fluid_mass_flow_rate / density
   cross_sectional_area = pipe_area_A_DL(None, pipe_inner_diameter, 1)
   return  volume / cross_sectional_area
