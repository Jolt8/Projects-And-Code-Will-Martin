from CoolProp.CoolProp import PropsSI, Props
import CoolProp.CoolProp as cp

import colebrook

import numpy as np

import cantera as ct
from cantera import ThermoPhase

import math 

import scipy
from scipy import optimize

import fluids 
from fluids import two_phase

from pyfluids import Fluid, FluidsList, Input
import scipy.optimize 



"""

state_values = [self.temperature, self.pressure, self.density, self.enthalpy, self.entropy, self.internal_Energy]
state_letters = ["T", "P", "D", "H", "S", "U"]

non_empty_props = {k: v for k, v in state_values if v is not None}

"""

def test2():
    temperature = 300
    pressure = 100000
    density = None
    enthalpy = None
    entropy = None
    internal_energy = None
    X = "CO2:1"
    state_values = [temperature, density, pressure, enthalpy, entropy, internal_energy]
    state_letters = ["T", "D", "P", "H", "S", "U"]

    non_empty_props = []
    for i in range(len(state_values)): 
        if state_values[i] is not None:
            non_empty_props.append(i)
        else:
            pass
    
    fluid = ct.Solution("gri30.yaml")
    
    prop1 = state_letters[non_empty_props[0]]
    prop2 = state_letters[non_empty_props[1]]
    val1 = state_values[non_empty_props[0]]
    val2 = state_values[non_empty_props[1]]
    method_call = f"fluid.{prop1}{prop2} = ({val1}, {val2})"
    
    exec(method_call)
    fluid.X = X
    
    print(fluid.HP)
    
def test3(input_file, X, temperature, pressure, density, enthalpy, entropy, internal_energy):
    fluid = ct.Solution(input_file)
    fluid.X = X
    
    if temperature is not None and pressure is not None:
        fluid.TP = temperature, pressure
    elif temperature is not None and density is not None:
        fluid.TD = temperature, density
    elif pressure is not None and enthalpy is not None:
        fluid.PH = pressure, enthalpy
    elif pressure is not None and entropy is not None:
        fluid.PS = pressure, entropy
    elif pressure is not None and internal_energy is not None:
        fluid.PU = pressure, internal_energy




class stream:
    """
    stream
        information about the working fluid in each point
            fluid
            quality
            temperature
            pressure
            enthalpy
            entropy
            gibbs free energy
            heat of combustion
            mass flow
            etc
    """
   
    def __init__(self, fluid_name, temperature, pressure, quality, mass_flow):
        self.fluid_name = fluid_name
        self.fluid.T = temperature
        self.pressure = pressure
        self.density = density 
        self.enthalpy
        self.entropy 
        
        self.fluid = ct.solution(self.fluid_name)
        
        
        self.mass_flow = mass_flow
        
    def set_state(self):
        
        if self.temperature is not None:
            self.fluid.T = self.temperature
        if self.pressure is not None:
            self.fluid.P = self.pressure
        if self.density is not None:
            self.fluid.D = self.density
        if self.enthalpy is not None:
            self.fluid.H = self.enthalpy
        if self.entropy is not None:
            self.fluid.S = self.entropy
        if self.internal_energy is not None:
            self.fluid.U = self.internal_energy
        if self.X is not None:
            self.fluid.U = self.internal_energy
        
    
        
        
            
        
        
    
    
        
    
    def calculations(self):
        self.enthalpy = cp.PropsSI("H", "T", self.temperature, "P", self.pressure, self.fluid)
        self.entropy = cp.PropsSI("S", "T", self.temperature, "P", self.pressure, self.fluid)
        self.gibbs_free_energy = self.enthalpy + self.entropy
        self.quality = find_Q_final_with_Tsat(saturation_temperature=self.temperature, actual_pressure=self.pressure)





state_values = [self.temperature, self.pressure, self.density, self.enthalpy, self.entropy, self.internal_Energy]
state_letters = ["T", "P", "D", "H", "S", "U"]

non_empty_props = {k: v for k, v in state_values if v is not None}







    
    
    




