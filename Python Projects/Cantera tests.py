from CoolProp.CoolProp import PropsSI, Props
import CoolProp.CoolProp as cp

import colebrook

import numpy as np

import cantera as ct
from cantera import ThermoPhase

import math 

import scipy
from scipy import optimize

import scipy.optimize 

#gas1 = ct.Solution("gri30.yaml")


fluid = ct.Solution("air.yaml")

# Set a temperature (in K) and find the saturation pressure
temperature = 450  # Boiling point of water at 100°C
fluid.TPX = temperature, cp.PropsSI("P", "T", temperature, "Q", 0, "Water"), "O2:1" # Set temperature and vapor quality (0 = saturated liquid)

# Get the saturation pressure (vapor pressure at this temperature)
vapor_pressure = fluid.P

test_fluid1 = ct.Solution("air.yaml")

print(test_fluid1.species)

test_fluid1.TPX = (300, 300000, "O2:1")


print(test_fluid1.DPX)
#test_fluid1.Q = 0.5
#print(test_fluid1.Q)
#ThermoPhase()
"""
fluid2 = ct.Solution("liquidvapor.yaml", "liquid-water-IAPWS95")

p_sat = fluid2.P_sat
#fluid2.TPX = (300, 300000, "H2O:1")
fluid2.TP = (300, 100000)
enthalpy = fluid2.h
print(enthalpy)
fluid2.h = enthalpy + 100000
fluid2.q = (1)
fluid2.phase_of_matter = ("gas")
#fluid2.v = 1


fluid2()
print(fluid2)
"""

fluid2 = ct.Solution("gri30.yaml")

fluid2.TPX = (300, 100000, "H2O:0.5 CH4:0.5")
fluid2.equilibrate("TP")

print(fluid2.TP)
print(fluid2.Y)
print(fluid2.density)


fluid3 = ct.Water()

fluid3.TPX = (300, 100000, "H2O:1")
print(fluid3.h)


fluid3.TPX = (300, 1000000, "H2O:1")
print(fluid3.h)

"""
print(f"Vapor pressure of water at {temperature} K: {vapor_pressure / 1e5:.2f} bar")
print(cp.PropsSI("P", "T", temperature, "Q", 0, "Water"))

print(cp.PropsSI("P", "T", 273 + 20, "Q", 0, "Butane"))


water = ct.Water()


fluid1 = ct.Solution("gri30.yaml")
fluid1.TPX = (300, 100000, "CH4:1")
q1 = ct.Quantity(fluid1, mass=5)

fluid2 = ct.Solution("gri30.yaml")
fluid2.TPX = (400, 200000, "O2:1")
q2 = ct.Quantity(fluid2, mass=5)


fluid3 = q1 + q2
print(fluid3.TPY)

# Given temperature and pressure
temperature = 300  # K
pressure = 90000    # 5 bar

# Set state with temperature and pressure
water.TP = temperature, pressure
print(water.Q)

if water.phase_of_matter == "liquid":
    quality = 0.0
elif water.phase_of_matter == "gas":
    quality = 1.0
else:


h = water.h  
water.TQ = temperature, 0
h_liquid = water.h
water.TQ = temperature, 1
h_vapor = water.h

quality = (h - h_liquid) / (h_vapor - h_liquid)

print(f"Quality: {quality:.4f}")
"""