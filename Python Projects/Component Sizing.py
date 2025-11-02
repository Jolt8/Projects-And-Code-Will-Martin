from CoolProp.CoolProp import PropsSI, Props
import CoolProp.CoolProp as cp

import numpy as np

import math 

import cantera as ct

import scipy 
from scipy import optimize

from pyfluids import Fluid, FluidsList, Input

import ht as ht

import networkx as nx

fluid = "Isobutane"
heat_wattage = 400
temp1 = 300
temp2 = 400
#enthalpy1 = cp.PropsSI("H", "T", temp1, "Q", 0, fluid)
#enthalpy2 = cp.PropsSI("H", "T", temp2, "Q", 1, fluid)
enthalpy1 = cp.PropsSI("H", "T", temp1, "P", 500000, fluid)
enthalpy2 = cp.PropsSI("H", "T", temp2, "P", 1500000, fluid)

delta_H = enthalpy2 - enthalpy1 

mass_flow = heat_wattage / delta_H
print("mass flow", mass_flow)

#liquid_volumetric_flow = mass_flow / cp.PropsSI("D", "T", temp1, "Q", 0, fluid) 
#vapor__volumetric_flow = mass_flow / cp.PropsSI("D", "T", temp2, "Q", 1, fluid)
print("m3/kg temp1", cp.PropsSI("V", "T", temp1, "Q", 0, fluid))
print("m3/kg temp2",cp.PropsSI("V", "T", temp2, "Q", 1, fluid))
print("density temp1", cp.PropsSI("D", "T", temp1, "Q", 0, fluid))
print("density temp2", cp.PropsSI("D", "T", temp2, "Q", 1, fluid))

liquid_volumetric_flow = mass_flow / cp.PropsSI("D", "T", temp1, "P", 500000, fluid)
vapor__volumetric_flow = mass_flow / cp.PropsSI("D", "T", temp2, "P", 1500000, fluid)
print(cp.PropsSI("V", "T", temp1, "P", 500000, fluid), cp.PropsSI("V", "T", temp2, "P", 1500000, fluid))
print(cp.PropsSI("D", "T", temp1, "P", 500000, fluid), cp.PropsSI("D", "T", temp2, "P", 1500000, fluid))

print("liquid volumetric flow", liquid_volumetric_flow*1000, "vapour volumetric flow", vapor__volumetric_flow*1000)

def gear_pump(number_of_teeth, gear_width_to_outer_diameter_ratio, liters_min, rpm):
    cc_min = liters_min * 1000
    result_container = {}
    def objective(module_array):
        module = module_array[0]
        pitch_diameter = module * number_of_teeth
        outer_diameter = module * (2 + number_of_teeth)
        gear_width = 1
        #gear_width = gear_width_to_outer_diameter_ratio * outer_diameter
        result_container["gear_width"] = gear_width
        print(
            module,
            gear_width,
            outer_diameter,
            ((math.pi / 2) * gear_width * (outer_diameter**2 - pitch_diameter**2)) - (cc_min / rpm)
        )
        return (math.pi / 2) * gear_width * (outer_diameter**2 - pitch_diameter**2) - (cc_min / rpm)
    result = scipy.optimize.root(objective, [1.5])
    final_gear_width = result_container.get('gear_width', None)
    return result, final_gear_width


print(gear_pump(
    number_of_teeth=12,
    gear_width_to_outer_diameter_ratio=0.5,
    liters_min=0.25,
    rpm=100
))


#def test2(*args):
    #i = args.index(None)
    #return scipy.optimize.fsolve(lambda x: gear_pump(*args[:i], x, *args[i+1:]), 0.5)
#print("Final test", test2(12, 0.5, 0.25, None))


def vane_expander(liters_min, rpm, rotor_width_to_diameter_ratio, rotor_diameter_to_stator_diameter_ratio):
    def objective(rotor_diameter_array):
        rotor_diameter = rotor_diameter_array[0]
        rotor_width = rotor_diameter* rotor_width_to_diameter_ratio 
        






