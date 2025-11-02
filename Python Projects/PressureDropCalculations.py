from CoolProp.CoolProp import PropsSI, Props
#from CoolProp.CoolProp import IProps, get_Fluid_index
#from CoolProp.CoolProp import PhaseSI
import CoolProp.CoolProp as cp

import colebrook

import numpy as np

import sys

import math 

import scipy
from scipy import optimize

import fluids 
from fluids import two_phase

from pyfluids import Fluid, FluidsList, Input
import scipy.optimize 

import BasicThermalCalculations
from BasicThermalCalculations import find_Q_final_with_Tsat

import cantera as ct

#The dynamic_viscosity of water is around 0.00089
def reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, characteristic_length):
    """Finds the Reynold's number of a fluid in a pipe 
        The Reynolds number predicts whether or not flow in a pipe will be laminar or not
        Re < 2300 = turbulent
        Re > 2300 = Laminar

    Args:
        Re (dimensionless): Reynold's number
        fluid (str): 
        temp (K): temperature of the fluid
        pressure (p): pressure inside the pipes
        fluid_density (kg/m3): _description_
        fluid_velocity (m/s): 
        diameter (m): inner diameter of the pipe

    Returns:
        Any variable when substituted with "none"
        can't solve for fluids
    """
    fluid_density = cp.PropsSI("D", "T", fluid_temp, "P", fluid_pressure, fluid)
    dynamic_viscosity = cp.PropsSI("V", "T", fluid_temp, "P", fluid_pressure, fluid)
    return (fluid_density * fluid_velocity * characteristic_length) / (dynamic_viscosity)
#print ("reynold's number", reynolds_number("Water", 300, 100000, 2, 0.1))

def cantera_reynolds_number(solution, fluid_velocity, inner_characteristic_length):
    """Finds the Reynold's number of a fluid in a pipe 
        The Reynolds number predicts whether or not flow in a pipe will be laminar or not
        Re < 2300 = turbulent
        Re > 2300 = Laminar

    Args:
        solution (obj): cantera solution class
        Re (dimensionless): Reynold's number
        inner_characteristic_length (m): inner diameter of the pipe
    """
    return (solution.density * fluid_velocity * inner_characteristic_length) / (solution.viscosity)
#print ("reynold's number", reynolds_number("Water", 300, 100000, 2, 0.1))


def pipe_pressure_drop(fluid, fluid_temp, fluid_pressure, mass_flow, pipe_diameter, pipe_length, absolute_pipe_roughness):
    """Finds Pressure loss with Darcy Weisback Friction Factor

    Args:
        fluid (str): Coolprop compatible fluid
        fluid_temp (K): temperature of the fluid in the pipe
        fluid_pressure (p): pressure of the fluid in the pipe to find the density 
        fluid_velocity (m/s): 
        pipe_diameter (m): 
        pipe_length (m): 
        pipe_roughness (mm): Absolute roughness of the pipe (Copper pipe is around 0.0015mm)
            Note that the colebrook module uses relative roughness which is the absolute roughness divided by the pipe's diameter

    Returns:
        Pressure drop (p): total pressure drop 
    """
    pipe_area = (math.pi * (pipe_diameter / 2) ** 2)
    fluid_density = cp.PropsSI("D", "T", fluid_temp, "P", fluid_pressure, fluid)
    volumetric_flow = mass_flow / fluid_density
    fluid_velocity = volumetric_flow / pipe_area
    Re = reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_diameter)
    relative_pipe_roughness = (absolute_pipe_roughness / pipe_diameter)
    friction_factor = colebrook.bntFriction(Re, relative_pipe_roughness, 4)
    #print (friction_factor)
    #print (fluid_density)
    return friction_factor * (pipe_length / pipe_diameter) * ((fluid_density * fluid_velocity ** 2) / 2)

#print ("Pipe_pressure_drop", pipe_pressure_drop("Water", 300, 100000, 1, 0.1, 1, 0.0015))


def cantera_pipe_pressure_drop(solution, mass_flow, inner_characteristic_length, pipe_length, absolute_pipe_roughness):
    """Finds Pressure loss with Darcy Weisback Friction Factor
    REMINDER: absolute pipe roughness is different from relative pipe roughness

    Args:
        fluid (str): Coolprop compatible fluid
        fluid_temp (K): temperature of the fluid in the pipe
        fluid_pressure (p): pressure of the fluid in the pipe to find the density 
        fluid_velocity (m/s): 
        pipe_diameter (m): 
        pipe_length (m): 
        pipe_roughness (mm): Absolute roughness of the pipe (Copper pipe is around 0.0015mm)
            Note that the colebrook module uses relative roughness which is the absolute roughness divided by the pipe's diameter

    Returns:
        Pressure drop (p): total pressure drop 
    """
    gravitational_acceleration = 9.8
    volumetric_flow = mass_flow / solution.density
    
    pipe_area = (math.pi * inner_characteristic_length * pipe_length)
    fluid_velocity = volumetric_flow / pipe_area
    
    Re = cantera_reynolds_number(solution, fluid_velocity, inner_characteristic_length)
    relative_pipe_roughness = (absolute_pipe_roughness / inner_characteristic_length)
    friction_factor = colebrook.bntFriction(Re, relative_pipe_roughness, 4)
    #print (friction_factor)
    #print (fluid_density)
    return friction_factor * (pipe_length / inner_characteristic_length) * (solution.density / 2) * fluid_velocity

#print ("cantera Pipe_pressure_drop", cantera_pipe_pressure_drop("gri30.yaml", "H2O:1", 300, 100000, 1, 0.1, 1, 0.0015))

def head_loss(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_diameter, pipe_length, pipe_roughness, gravitational_acceleration):
    Re = reynolds_number(fluid, fluid_temp, fluid_pressure, fluid_velocity, pipe_diameter)
    friction_factor = colebrook.bntFriction(Re, pipe_roughness, 4)
    return friction_factor * (pipe_length/pipe_diameter) * ((fluid_velocity ** 2) / (2 * gravitational_acceleration))


def cantera_head_loss(solution, fluid_velocity, inner_characteristic_length, pipe_length, pipe_roughness):
    gravitational_acceleration = 9.8
    Re = cantera_reynolds_number(solution, fluid_velocity, inner_characteristic_length)
    friction_factor = colebrook.bntFriction(Re, pipe_roughness, 4)
    return friction_factor * (1 / (2 * gravitational_acceleration)) * (pipe_length / inner_characteristic_length) * (fluid_velocity ** 2)

#print ("Head loss", head_loss("Water", 300, 100000, 1, 0.1, 1, 0.0015, 9.8))

#print(two_phase.two_phase_dP(m=0.6, x=0.1, rhol=915., rhog=2.67, mul=180E-6, mug=14E-6,
#sigma=0.0487, D=0.05, L=1.0))

#print(cp.PropsSI("T", "P", 300000, "Q", 0.1, "Isobutane"))
#print(find_Q_final_with_Tsat("Isobutane", 300, 300000))

def incremental_pressure_drop(fluid, fluid_in_temp, fluid_in_pressure, fluid_mass_flow, pipe_diameter, pipe_length, absolute_pipe_roughness):
    iterations = 10
    changing_fluid_pressure = fluid_in_pressure
    total_pressure_drop = 0

    for i in range(iterations):
        # Recalculate segment length if pipe_length has been adjusted
        pipe_segment_length = pipe_length / iterations
        
        quality = find_Q_final_with_Tsat(fluid, fluid_in_temp, changing_fluid_pressure)
        print(quality)
        liquid_density = cp.PropsSI("D", "T|liquid", fluid_in_temp, "P", changing_fluid_pressure, fluid)
        
        try:
            gas_density = cp.PropsSI("D", "T|gas", fluid_in_temp, "P", changing_fluid_pressure, fluid)
        except:
            gas_density = None
        
        liquid_viscosity = cp.PropsSI("V", "T|liquid", fluid_in_temp, "P", changing_fluid_pressure, fluid)
        
        try:
            gas_viscosity = cp.PropsSI("V", "T|gas", fluid_in_temp, "P", changing_fluid_pressure, fluid)
        except:
            gas_viscosity = None
        
        critical_pressure = cp.PropsSI("p_critical", fluid)
        #surface_tension = cp.PropsSI("I", "P", changing_fluid_pressure, "Q", quality, fluid)
        
        
        #print(fluid_mass_flow, quality, fluid_in_pressure, changing_fluid_pressure, liquid_density, pipe_diameter, pipe_segment_length, 
            #gas_density, liquid_viscosity, gas_viscosity, surface_tension, 
            #critical_pressure, absolute_pipe_roughness)
        
        #if quality == 0 or quality == 1:
        pressure_drop_for_segment = two_phase.two_phase_dP(
            fluid_mass_flow, quality, liquid_density, pipe_diameter, pipe_segment_length, 
            gas_density, liquid_viscosity, gas_viscosity, sigma=None, 
            P=changing_fluid_pressure, Pc=critical_pressure, roughness=absolute_pipe_roughness, Method=None,
        )
        print(total_pressure_drop, pressure_drop_for_segment, changing_fluid_pressure)
        """else:
            pressure_drop_for_segment = two_phase.Kim_Mudawar(
                fluid_mass_flow, quality, liquid_density, gas_density, liquid_viscosity, gas_viscosity, surface_tension, pipe_diameter, pipe_length
            )"""
        
        #print(pressure_drop_for_segment, changing_fluid_pressure, pipe_length, pipe_segment_length, quality, i)
        
        changing_fluid_pressure = changing_fluid_pressure - pressure_drop_for_segment
        total_pressure_drop = total_pressure_drop + pressure_drop_for_segment
        
        # If pressure goes below threshold, adjust pipe_length (and thus segment length) for subsequent iterations
        """
        if changing_fluid_pressure < 100000:
            changing_fluid_pressure = 100000
            pipe_length /= 10
            pipe_segment_length /= 10
        elif changing_fluid_pressure > 3000000:
            changing_fluid_pressure = 100000
            pipe_length *= 10
            pipe_segment_length *= 10
        """
            
    return total_pressure_drop

#print("incremental_pressure_drop", incremental_pressure_drop("Isobutane", 300, 300000, 0.30, 0.05, 0.1, 0.00015))


def cantera_incremental_pressure_drop(solution, mass_flow, pipe_diameter, pipe_length, absolute_pipe_roughness):
    iterated_solution = ct.Water()
    iterated_solution.TPX = solution.TPX
    
    iterations = 10
    total_pressure_drop = 0
    
    for i in range(iterations):
        # Recalculate segment length if pipe_length has been adjusted
        pipe_segment_length = pipe_length / iterations
        
        liquid_viscosity = solution.viscosity
        
        critical_pressure = solution.critical_pressure
        #print(iterated_solution.P)
        pressure_drop_for_segment = two_phase.two_phase_dP(
            mass_flow, solution.Q, solution.density, pipe_diameter, pipe_segment_length, 
            iterated_solution.density, liquid_viscosity, iterated_solution.viscosity, sigma=None, 
            P=iterated_solution.P, Pc=iterated_solution.critical_pressure, roughness=absolute_pipe_roughness, Method=None,
        )
        iterated_solution.TP = None, (iterated_solution.P - pressure_drop_for_segment)
        total_pressure_drop = total_pressure_drop + pressure_drop_for_segment
        
        
        # If pressure goes below threshold, adjust pipe_length (and thus segment length) for subsequent iterations
        """
        if changing_fluid_pressure < 100000:
            changing_fluid_pressure = 100000
            pipe_length /= 10
            pipe_segment_length /= 10
        elif changing_fluid_pressure > 3000000:
            changing_fluid_pressure = 100000
            pipe_length *= 10
            pipe_segment_length *= 10
        """
            
    return total_pressure_drop

solution = ct.Water()
solution.TP = 300, 500000
#print("cantera_incremental_pressure_drop", cantera_incremental_pressure_drop(solution, 0.05, 0.01, 100, 0.00015))



def cantera_optimize_pipe_length(target_pressure_drop, pipe_length_guess, solution, mass_flow, pipe_diameter, absolute_pipe_roughness):
    def objective(pipe_length):
        result = cantera_incremental_pressure_drop(
        solution, mass_flow, pipe_diameter, pipe_length[0], absolute_pipe_roughness
        )
        #print(f"Objective evaluation for pipe length {pipe_length}: Pressure drop result: {result}")
        return abs(result - target_pressure_drop)
    optimized = None
    
    retry_count = 0
    max_retries = 500  # Prevent infinite loops by limiting retries

    while optimized is None and retry_count < max_retries:
        try:
            optimized = scipy.optimize.minimize(objective, [pipe_length_guess], method='Nelder-Mead')
        except:
            print("GUESS TO HIGH")
            pipe_length_guess /= 10  # Reduce the guess if an error occurs
            retry_count += 1  # Track how many retries have been done
    
    if retry_count >= max_retries:
        print("Max retries reached. Optimization did not succeed.")
    
    return optimized
"""
test_solution = ct.Water()
test_solution.TPX = 300, ct.one_atm, "H2O:1"
optimized_length = cantera_optimize_pipe_length(
    target_pressure_drop=30000,         
    pipe_length_guess=1,                  
    solution=test_solution,
    mass_flow=0.05,                      # 0.05 kg/s — higher flow rate
    pipe_diameter=0.01,                 # 0.02 m — 20 mm inner diameter
    absolute_pipe_roughness=0.00015      # 1.5 µm — smooth copper pipe
)

print(optimized_length)
"""

"""
def optimize_pipe_length(target_pressure_drop, initial_guess, fluid, fluid_in_temp, fluid_in_pressure, fluid_mass_flow, pipe_diameter, absolute_pipe_roughness):
    def objective(pipe_length):
        result = incremental_pressure_drop(
        fluid, fluid_in_temp, fluid_in_pressure, fluid_mass_flow, pipe_diameter, pipe_length[0], absolute_pipe_roughness
        )
        #print(f"Objective evaluation for pipe length {pipe_length}: Pressure drop result: {result}")
        return abs(result - target_pressure_drop)
    optimized = None
    
    retry_count = 0
    max_retries = 500  # Prevent infinite loops by limiting retries

    while optimized is None and retry_count < max_retries:
        try:
            optimized = scipy.optimize.minimize(objective, [initial_guess], method='Nelder-Mead')
        except:
            initial_guess /= 10  # Reduce the guess if an error occurs
            retry_count += 1  # Track how many retries have been done
    
    if retry_count >= max_retries:
        raise ValueError("Max retries reached. Optimization did not succeed.")
    
    return optimized


optimized_length = optimize_pipe_length(
    target_pressure_drop=75000,         
    initial_guess=1.0,                  
    fluid='Isobutane',                  
    fluid_in_temp=350,                  # 350 K — typical superheat temp
    fluid_in_pressure=800000,           # 800 kPa — high side ORC pressure
    fluid_mass_flow=0.05,               # 0.05 kg/s — higher flow rate
    pipe_diameter=0.01,                 # 0.02 m — 20 mm inner diameter
    absolute_pipe_roughness=1.5e-6      # 1.5 µm — smooth copper pipe
)

print("non-cantera version", optimized_length)
"""



