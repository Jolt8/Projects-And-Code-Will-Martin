from CoolProp.CoolProp import PropsSI, Props
import CoolProp.CoolProp as cp

import numpy as np

import math 

import cantera as ct

import scipy 
from scipy import optimize

from pyfluids import Fluid, FluidsList, Input
import scipy.optimize 

import ht as ht

from PressureDropCalculations import pipe_pressure_drop, incremental_pressure_drop, cantera_pipe_pressure_drop, cantera_incremental_pressure_drop

from HeatExchangers import cantera_air_heat_loss, heat_duty, cantera_overall_heat_transfer_coeff, heat_transfer_coeff, cantera_ht_boiling_flow_chen_bennet, cantera_ht_external_pipe, cantera_ht_external_plate, cantera_ht_internal_conv, cantera_ht_submerged_coil

import networkx as nx

"""
coolprop functions to know
    
    non trivial
        Dmolar
        Hmolar
        Smolar
        molarmass
        
        
    trivial
        molar mass - cp.PropsSI("molarmass", "fluid")
    
Note:
    internal energy - energy from molecular motion, latent heat, and chemical energy
    enthalpy - just from heat energy and the substance's state of matter
"""


"""
inputs that will be required coma chemical reaction
    reaction rate
    input concentrations
    output concentration 
    products
    reactants 
    heat of reaction
    equilibrium constant (if it's an equilibrium reaction)
    mass flow rate of each reactant, mass flow rate of reactants combined (with an ideal or non-ideal ratio), product mass flow rate 
    desired concentration of product
"""

"""
additional chemical engineering functions that will be required:
    distilation towers (rigorous method)
    compressors, expanders, and pumps
    heat exchangers (already covered by me)
    different types of reactants
"""
"""
stuff I'm not really sure about
    is there an easy way to find the impact that temperature, pressure, and concentrations will have on the reaction
    should I use classes?
    should I use a seperate class for equilibrium and non-equilibrium reactions
    should I use a seperate class for reactants and products?


"""
"""
classes needed
    compnent
        calculate (inlet_temp, pressure, mass flow --> outlet temp, pressure, mass flow)
        energy balance
        pressure drop
    thermodynamic components
        things that apply to all
            efficiency
            mass flow
            diameter
            heat exchange with environment
            roughness and pressure drop
        
        
        compressor 
            efficiency
        expander
            efficiency
        pump
            efficiency
            mass flow
            
        pipe
            mass flow
            diameter
            pressure drop
            heat exchange with environment
            
        heatexchanger
        condenser
        evaporator
    chemical engineering
        chemical reactor
        distillation column
        mixer
        splitter
        seperator
    supporting classes
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
                etc.
            there must be a way to split the stream class so it can handle splits and merges
            
        reaction
    system level classes
        system 
            takes list of component objects and initial conditions
            methods:
                run
                convergence_check




"""



#utility functions

#fuck this stupid class
""" 
class StreamManager:
    #should I just make the create_stream and other functions global utilities and have a seperate mass flow rate / transport class?
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(StreamManager, cls).__new__(cls)
            cls._instance.streams = {}
            cls._instance.mass_flow = {}
            cls._instance.volumetric_flows = {}
            cls._instance.molar_flows = {}
        return cls._instance

    def create_stream(self, file, cantera_name, stream_name, mass_flow=0, volumetric_flow=0, molar_flow=0):
        if stream_name in self.streams:
            raise ValueError(f"Stream '{stream_name}' already exists.")
        self.streams[stream_name] = ct.Solution(file, cantera_name)
        density = (self.streams[stream_name]).density
        molar_density = (self.streams[stream_name]).density_mole
        
        if (mass_flow != 0) or (mass_flow is not None):
            #for future me: the ^ operator just means xor
            self.mass_flow[stream_name] = mass_flow
            self.volumetric_flows[stream_name] = mass_flow / density
            self.molar_flows[stream_name] = mass_flow / molar_density
        elif (volumetric_flow != 0) or (volumetric_flow is not None):
            self.volumetric_flows[stream_name] = volumetric_flow
            self.mass_flow[stream_name] = volumetric_flow * density
            self.molar_flows[stream_name] = self.mass_flow[stream_name] / molar_density
        elif (molar_flow != 0) or (molar_flow is not None):
            self.molar_flows[stream_name] = molar_flow
            self.mass_flow[stream_name] = molar_flow * molar_density
            self.volumetric_flows[stream_name] = self.mass_flow[stream_name] / density
        elif mass_flow == 0 and volumetric_flow == 0 and molar_flow == 0:
            raise ValueError(f"At least one flow rate must be provided for '{stream_name}'")


    def get_stream(self, stream_name):
        if stream_name not in self.streams:
            raise ValueError(f"Stream '{stream_name}' not found.")
        return self.streams[stream_name]
    
    def show_streams(self):
        return self.streams
    
    def delete_stream(self, stream_name):
        if stream_name not in self.streams:
            raise ValueError(f"Stream '{stream_name}' not found.")
        del self.streams[stream_name]
        del self.mass_flow[stream_name]
        del self.volumetric_flows[stream_name]
        del self.molar_flows[stream_name]
    #please do note that splitting a stream is different from seperating a stream
        #splitting - splitting the mass flow of a stream
        #separating - separating the elements and compounds of a stream that also results in a split of mass flows between the two 
    def split_stream(self, stream_name, split_ratios):
        if stream_name not in self.streams:
            raise ValueError(f"Stream '{stream_name}' not found.")
        print("EeEE", split_ratios)
        try:
            total_ratio = sum(split_ratios.values())
        except:
            total_ratio = 1
        if not (0.99 <= total_ratio <= 1.01):
            raise ValueError("Split ratios must sum to 1.")
        original_mass_flow = self.mass_flow[stream_name]
        original_volumetric_flow = self.volumetric_flows[stream_name]
        original_molar_flow = self.molar_flows[stream_name]
        for new_name, ratio in split_ratios.items():
            if new_name in self.streams:
                raise ValueError(f"Stream '{new_name}' already exists.")
            self.streams[new_name] = self.streams[stream_name]
            self.mass_flow[new_name] = original_mass_flow * ratio
            self.volumetric_flows[new_name] = original_volumetric_flow * ratio
            self.molar_flows[new_name] = original_molar_flow * ratio
    def mix_streams(self, new_name, stream_names):
        if new_name in self.streams:
            raise ValueError(f"Stream '{new_name}' already exists.")
        
        total_mass_flow = sum(self.mass_flow[stream_name] for stream_name in stream_names)
        total_volumetric_flow = sum(self.volumetric_flows[stream_name] for stream_name in stream_names)
        total_molar_flow = sum(self.molar_flows[stream_name] for stream_name in stream_names)
        
        if total_mass_flow == 0:
            raise ValueError("Total mass flow is zero; cannot mix streams.")
        if total_volumetric_flow == 0:
            raise ValueError("Total volumetric flow is zero; cannot mix streams.")
        if total_molar_flow == 0:
            raise ValueError("Total molar flow is zero; cannot mix streams.")

        quantities = []
        for stream_name in stream_names:
            quantities.append(ct.Quantity(stream_names[stream_name], mass=(self.mass_flow[stream_name])))
        
        mixed_stream = sum(quantities)
        
        self.streams[new_name] = mixed_stream
"""
        
    
# Example usage
"""
if __name__ == "__main__":
    manager = StreamManager()

    manager.create_stream('air1', 'gri30.yaml', mass_flow=10)
    print(f"molar with mass: {manager.molar_flows['air1']} m3/s")
    manager.create_stream('air2', 'gri30.yaml', 0, 0, molar_flow=30)
    print(f"volumetric with moles: {manager.volumetric_flows['air2']} m3/s")
    
    
    #manager.split_stream('air', {'air1': 0.6, 'air2': 0.4})

    manager.create_stream('fuel', 'gri30.yaml', mass_flow=5)
    #manager.create_stream('air', 'gri30.yaml', mass_flow=10)
    manager.mix_streams('combustion_mix', ['air1', 'fuel'])

    mix = manager.get_stream('combustion_mix')
    print(mix.TPX)
    print(f"Mass flow of combustion mix: {manager.mass_flow['combustion_mix']} kg/s")
    print(f"Mass flow of combustion mix: {manager.volumetric_flows['combustion_mix']} m3/s")
"""
"""
component inputs (for all components)
    downstream_components
    component_name
    stream_name
    

"""

class Walls():
    def __init__(self, connected_components:list, ):
        #make sure that connect is in the right order if using a fluid to transport heat 
        self.connected_components = connected_components
        
        
class ConductanceHeatMover():
    def __init__(self, connected_components:list):
        Walls.__init__(self, connected_components)
        

class FluidHeatMover():
    def __init__(self, connected_components:list):
        Walls.__init__(self, connected_components)
        

class OpenHeatMover():
    def __init__(self, connected_components:list):
        Walls.__init__(self, connected_components)
        
        
class OpenHeatMover():
    def __init__(self, connected_components:list):
        Walls.__init__(self, connected_components)
        


def component_mix_streams(input_components:list, output_component):
    """
        example input:
        split_streams(splitter1, {compressor3: 0.5, compressor9: 0.5})
        split_streams(pipe1, {pipe2: 0.25, pipe3: 0.75})
        
        here's the referenced cantera docs:
        https://cantera.org/3.1/examples/python/reactors/mix1.html 
    """
    resevoirs = []
    
    mixer = ct.IdealGasReactor(None, name="Mixer")
    
    total_mass_flow = sum(input_components[i].mass_flow for i in range(len(input_components)))
    
    for i in range(len(input_components)):
        resevoirs.append(ct.resevoir(input_components[i].streams, str("resevoir ", i)))
        ct.MassFlowController(resevoirs[i], mixer, m_dot=input_components[i].mass_flows, name = str("controller ", i))
        
    output_resevoir = ct.resevoir(None, None)
    
    outlet = ct.Valve(mixer, output_resevoir, total_mass_flow)
    
    sim = ct.ReactorNet([mixer])
    
    sim.advance_to_steady_state()
    compressor.__init__
    output_component.inlet_fluids[0].TPX = mixer.Thermo
    output_component.mass_flow = mixer.mass_flow_rates
    

"""
print("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEE", input_component.inlet_fluids)
print(input_component.inlet_fluids[0].TPX)
print(input_component.inlet_fluids[1].TPX)


A = ct.Quantity(input_component.inlet_fluids[0], constant="HP")
quantities = ct.Quantity((input_component.input_components[i] for i in range(len(input_component.inlet_fluids))), constant="HP", )
input_component.inlet_fluids = input_component.inlet_fluids[0] + input_component.inlet_fluids[1]
input_component.mass_flow = mixer.mass_flow_rates
"""


    

    



class component_template():
    def __init__(self, upstreams:list, name, inlet_fluids:list, mass_flow, inner_characteristic_length, outer_characteristic_length):
        #how do I make it so that this component doesn't require stream inputs from the user but 
        self.upstreams = upstreams
        self.name = name
        self.inlet_fluids = inlet_fluids
        self.mass_flow = mass_flow
        self.inner_characteristic_length = inner_characteristic_length
        self.outer_characteristic_length = outer_characteristic_length
        self.thickness = self.inner_characteristic_length - self.outer_characteristic_length
        
        try:
            for i in range(len(self.inlet_fluids)):
                density = inlet_fluids[i].density
                
                area = (math.pi * ((inner_characteristic_length / 2) ** 2))
                
                if self.mass_flow is not None:
                    self.volumetric_flow = mass_flow / density 
                    self.velocity = self.volumetric_flow / area
                elif self.volumetric_flow is not None:
                    self.mass_flow = self.volumetric_flow / density
                    self.velocity = self.volumetric_flow / area
                elif self.velocity is not None:
                    self.volumetric_flow = area / self.velocity
                    self.mass_flow = self.volumetric_flow / density
        except:
            pass
            
        #def environmental_heat_exchange():
            


class Expander(component_template):
    def __init__(self, upstreams:list, name, inlet_fluids:list, mass_flow, inner_characteristic_length, outer_characteristic_length, efficiency, pressure_change, enthalpy_change, work_in, compression_ratio):
        component_template.__init__(self, upstreams, name, inlet_fluids, mass_flow, inner_characteristic_length, outer_characteristic_length)
        self.efficiency = efficiency
        self.pressure_change = pressure_change
        self.enthalpy_change = enthalpy_change
        self.work_in = work_in
        self.compression_ratio = compression_ratio
        self.outlet_fluids = []
    
    def update_state(self):
        #outlet_fluid = ct.Solution(str(self.inlet_fluids.source), str(self.inlet_fluids.name))
        #if len(self.inlet_fluids) >= 2:
            #self.inlet_fluids = list[mix_streams(self)]
        #else:
            #pass
        
        outlet_fluid = ct.Solution("gri30.yaml", "gri30")
        outlet_fluid.TPX = self.inlet_fluids[0].TPX
        
        if self.pressure_change != 0: #^ self.pressure_change is not None:
            initial_temperature = outlet_fluid.T
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            outlet_fluid.SP = initial_entropy, (initial_pressure + self.pressure_change)
            new_enthalpy = outlet_fluid.h
            
            isentropic_enthalpy_change = new_enthalpy - initial_enthalpy
            self.enthalpy_change = isentropic_enthalpy_change / self.efficiency
            self.work_in = self.enthalpy_change * self.mass_flow
            outlet_fluid.HP = initial_enthalpy + self.enthalpy_change, (initial_pressure + self.pressure_change)
        #the enthalpy change one might be unnecessary 
        elif self.enthalpy_change != 0: #^ self.enthalpy_change is not None:
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            new_enthalpy = initial_enthalpy + self.enthalpy_change
            outlet_fluid.HP = new_enthalpy, None
            
            self.pressure_change = initial_pressure - outlet_fluid.P
            self.work_in = self.enthalpy_change * self.mass_flow
        elif self.work_in != 0: #^ self.work_in is not None:
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            self.enthalpy_change = self.work_in / self.mass_flow
            
            new_enthalpy = initial_enthalpy + self.enthalpy_change
            outlet_fluid.HP = new_enthalpy, None
            
            self.pressure_change = initial_pressure - outlet_fluid.P
        elif self.compression_ratio != 0: #^ self.compression_ratio is not None:
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            self.pressure_change = initial_pressure * self.compression_ratio - initial_pressure
            
            outlet_fluid.SP = initial_entropy, (initial_pressure + self.pressure_change)
            new_enthalpy = outlet_fluid.h
            
            isentropic_enthalpy_change = new_enthalpy - initial_enthalpy
            self.enthalpy_change = isentropic_enthalpy_change / self.efficiency
            self.work_in = self.enthalpy_change * self.mass_flow
            outlet_fluid.HP = initial_enthalpy + self.enthalpy_change, (initial_pressure + self.pressure_change)
        self.outlet_fluids.append(outlet_fluid)

class Pipe(component_template):
    def __init__(self, upstreams:list, name, inlet_fluids:list, mass_flow, inner_characteristic_length, outer_characteristic_length, 
                 length, do_heat_loss=False, pipe_thermal_conductivity=0, expected_temperature_rise=30, air_temperature_initial=300, air_temperature_final=320, 
                 inner_wall_temperature=320, outer_wall_temperature=280, 
                 do_pressure_drop=False, pressure_drop_method="simple", absolute_pipe_roughness=0, 
                 ):
        """For pressure drop method, your two options are "simple" and "twophase"
        
        
        """
        component_template.__init__(self, upstreams, name, inlet_fluids, mass_flow, inner_characteristic_length, outer_characteristic_length)
        self.outlet_fluids = []
        self.length = length
        
        self.do_heat_loss = do_heat_loss
        self.pipe_thermal_conductivity = pipe_thermal_conductivity
        self.expected_temperature_rise = expected_temperature_rise
        self.air_temperature_initial = air_temperature_initial
        self.air_temperature_final = air_temperature_final
        self.inner_wall_temperature = inner_wall_temperature
        self.outer_wall_temperature = outer_wall_temperature
        
        self.do_pressure_drop = do_pressure_drop
        self.absolute_pipe_roughness = absolute_pipe_roughness
        self.relative_pipe_roughness = absolute_pipe_roughness / self.inner_characteristic_length
        
        self.pipe_thermal_conductivity = pipe_thermal_conductivity
        
    def pressure_change(self):
        if self.pressure_drop_method == "simple":
            self.pressure_drop = cantera_pipe_pressure_drop(
                self.inlet_fluids[0], self.mass_flow, self.inner_characteristic_length, self.length, self.absolute_pipe_roughness)
            self.current_fluid.P = self.current_fluid.P - self.pressure_drop
        
        elif self.pressure_drop_method == "twophase":
            self.pressure_drop = cantera_incremental_pressure_drop(self.inlet_fluids[0], self.mass_flow, self.inner_characteristic_length, self.length, self.absolute_pipe_roughness)
            self.current_fluid.P = self.current_fluid.P - self.pressure_drop
            
    def heat_change(self):
        #FIXME: I need a good way to find the final temperatures of the internal fluid and the rise in air temperature
        self.heat_exchange = cantera_air_heat_loss(self.inlet_fluids[0], self.expected_temperature_rise, self.mass_flow, 
                                                   self.inner_wall_temperature, self.outer_wall_temperature, 
                                                   self.inner_characteristic_length, self.outer_characteristic_length, self.length,
                                                   self.absolute_pipe_roughness, self.pipe_thermal_conductivity, 
                                                   self.air_temperature_initial, self.air_temperature_final)
        time_spent_in_pipe = self.length / self.velocity 
        self.enthalpy_change = (self.heat_exchange * time_spent_in_pipe)
        self.current_fluid.H = self.currentfluid.H + (self.heat_exchange + self.current_fluid.H)
    
    def update_state(self):
        if self.do_pressure_drop == True:
            self.pressure_change()
        if self.do_heat_loss == True:
            self.heat_change()
            
            
class Splitter(component_template):
    """
    Example of outlet_stream_split variable: [0.8, 0.8, 9]
    """
    def __init__(self, upstreams:list, name, inlet_fluids:list, mass_flow, inner_characteristic_length, outer_characteristic_length, respective_split_ratios:list):
        component_template.__init__(self, upstreams, name, inlet_fluids, mass_flow, inner_characteristic_length, outer_characteristic_length)
        self.outlet_fluids = []
        
        

"""
    oh shit! how do we do heat / energy transfer?
        walls
            wall1(adiabadic_reactor, pipe1, wall_effectiveness)
                wall_effectiveness describes the ratio between the heat transfered into the other component and heat still lost to environment
                    ex. if a compressor or reactor was creating 500w of heat, a ratio of 0.5 would allow 250w to be actually transfered, the other 250w lost
                    if we didn't do this at all, we would lose 500w of heat regardless
        energy_transfer_fluid
            energy_transfer_fluid(adiabadic_reactor, pipe1, )
                for heat transfer using a working fluid
        wire / axle
            wire / axle(expander, pump)
"""


"""
start_state 
            start_state = fluid.TPX = ()
            start_position = (mixer, pipe1, (mass flow, molar flow, volumetric flow rate, or speed))
                component_after
            
        ambient_conditions
            air_temperature = ()
            air_pressure = ()
            
        chemical_reactions
            cantera_reactions = []
                #these will be used later in reactors
"""

class compressor(component_template):
    def __init__(self, upstreams:list, name, inlet_fluids:list, mass_flow, inner_characteristic_length, outer_characteristic_length, efficiency, pressure_change, enthalpy_change, work_in, compression_ratio):
        component_template.__init__(self, upstreams, name, inlet_fluids, mass_flow, inner_characteristic_length, outer_characteristic_length)
        self.outlet_fluids = []
        self.efficiency = efficiency
        self.pressure_change = pressure_change
        self.enthalpy_change = enthalpy_change
        self.work_in = work_in
        self.compression_ratio = compression_ratio
    
    def update_state(self):
        #outlet_fluid = ct.Solution(str(self.inlet_fluids.source), str(self.inlet_fluids.name))
        #if len(self.inlet_fluids) >= 2:
            #self.inlet_fluids = list[mix_streams(self)]
        #else:
            #pass
        
        outlet_fluid = ct.Solution(self.inlet_fluids[0].source)
        outlet_fluid.TPX = self.inlet_fluids[0].TPX
        
        if self.pressure_change != 0: #^ self.pressure_change is not None:
            initial_temperature = outlet_fluid.T
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            outlet_fluid.SP = initial_entropy, (initial_pressure + self.pressure_change)
            new_enthalpy = outlet_fluid.h
            
            isentropic_enthalpy_change = new_enthalpy - initial_enthalpy
            self.enthalpy_change = isentropic_enthalpy_change / self.efficiency
            self.work_in = self.enthalpy_change * self.mass_flow
            outlet_fluid.HP = initial_enthalpy + self.enthalpy_change, (initial_pressure + self.pressure_change)
        #the enthalpy change one might be unnecessary 
        elif self.enthalpy_change != 0: #^ self.enthalpy_change is not None:
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            new_enthalpy = initial_enthalpy + self.enthalpy_change
            outlet_fluid.HP = new_enthalpy, None
            
            self.pressure_change = initial_pressure - outlet_fluid.P
            self.work_in = self.enthalpy_change * self.mass_flow
        elif self.work_in != 0: #^ self.work_in is not None:
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            self.enthalpy_change = self.work_in / self.mass_flow
            
            new_enthalpy = initial_enthalpy + self.enthalpy_change
            outlet_fluid.HP = new_enthalpy, None
            
            self.pressure_change = initial_pressure - outlet_fluid.P
        elif self.compression_ratio != 0: #^ self.compression_ratio is not None:
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            self.pressure_change = initial_pressure * self.compression_ratio - initial_pressure
            
            outlet_fluid.SP = initial_entropy, (initial_pressure + self.pressure_change)
            new_enthalpy = outlet_fluid.h
            
            isentropic_enthalpy_change = new_enthalpy - initial_enthalpy
            self.enthalpy_change = isentropic_enthalpy_change / self.efficiency
            self.work_in = self.enthalpy_change * self.mass_flow
            outlet_fluid.HP = initial_enthalpy + self.enthalpy_change, (initial_pressure + self.pressure_change)
        self.outlet_fluids.append(outlet_fluid)

#water = ct.Solution("gri30.yaml", "gri30")
#water.TPX = 300, 100000, "H2O:1"


#compressor_1 = compressor([compressor_2], "compressor_1", water, 3, None, None, 0.80, 100000, 0, 0, 0, 0)

def mix_streams_with_upstreams(input_component):
    """
        example input:
        split_streams(splitter1, {compressor3: 0.5, compressor9: 0.5})
        split_streams(pipe1, {pipe2: 0.25, pipe3: 0.75})
        
        here's the referenced cantera docs:
        https://cantera.org/3.1/examples/python/reactors/mix1.html 
    """
    reservoirs = []
    
    total_mass_flow = sum(input_component.upstreams[i].mass_flow for i in range(len(input_component.upstreams)))

    for i in range(len(input_component.upstreams)):
        reservoirs.append(ct.Reservoir(input_component.upstreams[i].outlet_fluids[0], name=None))
    
    mixer = ct.IdealGasReactor(input_component.upstreams[0].outlet_fluids[0], name=None)
    
    for i in range(len(input_component.upstreams)):
        ct.MassFlowController(reservoirs[i], mixer, mdot=input_component.upstreams[i].mass_flow, name=None)
        
    output_reservoir = ct.Reservoir(input_component.upstreams[0].outlet_fluids[0])
    
    outlet = ct.Valve(mixer, output_reservoir, K=10)
    
    sim = ct.ReactorNet([mixer])
    
    sim.advance_to_steady_state()
    
    input_component.inlet_fluids = [mixer.thermo]
    input_component.mass_flow = total_mass_flow


class system():
    def __init__(self, starting_stream, starting_component, components):
        self.starting_stream = starting_stream
        self.starting_component = starting_component
        self.components = components
        self.total_mass_flow = self.starting_component.mass_flow
        #we should probably implement something like this
        #self.total_heat_wattage = total_heat_wattage
        
    
    """
    def process_comp_list(self):
        processed_comp_list = []
        starting_number = self.components.index(self.start_position)
        #starts at a certian position and loops over everything until it reaches the starting position again
        #ex. 3, 4, 5, 1, 2
        for i in range(starting_number, abs(starting_number - 2)):
            print(self.components[i].downstreams)
            print(len(self.components[i].downstreams))
            if len(self.components[i].downstreams) == 1:
                processed_comp_list.append(self.components[i])
            elif len(self.components[i].downstreams) > 1:
                processed_comp_list.append(self.components[i], "split", [self.components[splits] for splits in self.components[i].downstreams])
                i += len([self.components[i].downstreams])
            elif type(self.components).__name__ == splitter():
                processed_comp_list.append(self.components[i], "split", [self.components[splits] for splits in self.components[i].downstreams])
                i += len([self.components[i].downstreams])
            elif type(self.components).__name__ == distillation():
                processed_comp_list.append(self.components[i], "split", [self.components[splits] for splits in self.components[i].downstreams])
                i += len([self.components[i].downstreams])
            else:
                raise ValueError("downstreams input for ", self.component[i], "(", self.component[i].downstreams, ") is invalid")
            print("EEE", processed_comp_list)
        """
    
    def get_total_head_loss(self):
        total_head_loss = []
        for i in range(len(self.components)):
            if hasattr(self.components[i], "Pipe"):
                total_head_loss.append(self.components[i].head_loss)
    
    
    
    
    def split_streams(first_component, split_ratios:dict):
        """
            example input:
            split_streams(splitter1, {compressor3: 0.5, compressor9: 0.5})
            split_streams(pipe1, {pipe2: 0.25, pipe3: 0.75})
        """
        resevoirs = []
        i = 0
        for component, ratio in split_ratios:
            resevoirs.append(ct.resevoirs(split_ratios.items(component)), name=str(i))
            i += 1
        

                
    
    def better_run(self):
        self.get_total_head_loss = system.get_total_head_loss(self)
        #order the list starting from the starting position defined by the user
        self.components = self.components[self.components.index(self.starting_component):] + self.components[:self.components.index(self.starting_component)]
        
        G = nx.DiGraph()
        G.add_nodes_from(self.components)
        
        #you have to do this or else you will get an error that the dictionary changed during iteration 
        nodes = list(G.nodes)
        for node in nodes:
            G.add_edges_from([(node.upstreams[ds_index], node) for ds_index in range(len(node.upstreams))])
       
        visited = [] #change back to set() when done
        
        print([self.components[i].name for i in range(len(self.components))])
        self.starting_component.inlet_fluids = [self.starting_stream]
        self.starting_component.update_state() 
        visited.append(self.starting_component.name)
        
        print(self.starting_component.name)
        print(self.starting_component.inlet_fluids[0].T, self.starting_component.inlet_fluids[0].P, self.starting_component.mass_flow, self.starting_component.inlet_fluids[0].h)
        print(self.starting_component.outlet_fluids[0].T, self.starting_component.outlet_fluids[0].P, self.starting_component.mass_flow, self.starting_component.outlet_fluids[0].h)
        print("")
        
        
        
        start_index = self.components.index(self.starting_component)
        for src, node in nx.edge_bfs(G, self.components[start_index]):
            if node.name not in visited: #or len(node.upstreams) >= 2:
                src_count = sum(1 for s, _ in nx.edge_bfs(G, self.components[start_index]) if s == src)
                if nodes.index(node) == 0:
                    node.inlet_fluids = [self.starting_stream]
                    node.mass_flow = src.mass_flow
                #this is for mixing (the current component has two inputs)
                elif len(node.upstreams) >= 2:
                    mix_streams_with_upstreams(node)
                    #node.inlet_fluids = src.outlet_fluids
                    #node.mass_flow = sum(node.upstreams[i].mass_flow for i in range(len(node.upstreams)))
                #this is for splitting (the last component had multiple outputs)
                elif src_count >= 2:
                    node.inlet_fluids = src.outlet_fluids
                    node.mass_flow = src.mass_flow / src_count
                else:
                    node.inlet_fluids = src.outlet_fluids
                    node.mass_flow = src.mass_flow
                    
                
                """
                #for mixing
                elif len(node.inlet_fluids) >= 2:
                    node.inlet_fluids = src.inlet_fluids
                    node.mass_flow = sum(src.inlet_fluids) 
                #for 
                elif len(src.upstreams) >= 2:
                    node.inlet_fluids = src.outlet_fluids
                    node.mass_flow = src.mass_flow / len(src.upstreams)
                """
                
                    
                
                node.update_state()
                visited.append(node.name)
                
                print(node.name)
                print("current component previous connections:", list(node.upstreams[i].name for i in range(len(node.upstreams))))
                print(node.inlet_fluids[0].T, node.inlet_fluids[0].P, node.mass_flow, node.inlet_fluids[0].h, node.inlet_fluids)
                print(node.outlet_fluids[0].T, node.outlet_fluids[0].P, node.mass_flow, node.outlet_fluids[0].h, node.outlet_fluids)
                print(visited)
                print("")
                
        new_dict = []
        for src, node in nx.edge_bfs(G, self.components[start_index]):
            new_dict.append((src.name, "-->", node.name))
        print(new_dict)
        """
        #BFS to visit all nodes while handling splits
        for src, node in nx.bfs_edges(G, self.starting_component):
            if node not in visited:
                successors = dict(nx.bfs_successors(G, self.starting_component))
                predecessors = dict(nx.bfs_predecessors(G, self.starting_component))
                print("successors", successors)
                print("predecessors", predecessors)
                #print("edge test", nx.bfs_edges(G ))
                if len(successors[node]) >= 2:
                    print("hiiii", len(successors[node]), node.name)
                    try:
                        if node.split_ratios is not None and hasattr(node, "splitter"):
                            compressor
                            manager.split_stream(node.stream, ({node.stream + str(i) : split_ratio} for split_ratio, i in node.split_ratios))
                        #manager.split_stream('air', {'air1': 0.6, 'air2': 0.4}))
                        #manager.mix_streams('combustion_mix', ['air1', 'fuel'])
                    except:
                        manager.split_stream(node.stream, {node.stream + str(1) : 0.5, node.stream + str(2) : 0.5})
                        print(manager.streams)
                        print(manager.mass_flow["water"])
                        print(manager.mass_flow["water1"])
                        print(successors[node])
                        print(node.name)
                        manager.split_stream(node.stream, ({node.stream + str(i) : 0.5} for i in range(len(successors[node]))))
                                                   
                    
                #elif len([predecessors[node]]) >= 2:
                    #manager.mix_streams(node.name, predecessors[node])
                
                else:
                    print("hi")
                
                
                
                
                print(1)
                
                
                print("EEE", src, node)
                print("bfs_layers", dict(enumerate(nx.bfs_layers(G, self.start_position))))
                print("bfs_tree", list(nx.bfs_tree(G, self.start_position)))
                print("bfs_predecessors", dict(nx.bfs_predecessors(G, self.start_position)))
                print("bfs_successors", dict(nx.bfs_successors(G, self.start_position)))
                
                visited.add(node)
        """

"""
    I think the best way to structure this is something like this:
    
        whenever a split is encountered, create a new stream for this split and run it to the end where it either terminates on a product component (output)
        or a mix component
        then, do the next item in the split and do the same
        
        I think the way this will be structured is that we'll have each split denoted with x.x.x.x.x... structure that keeps track where it's positioned in the split network
        the only issue with that is that if a stream is briefly split and then mixed, should we keep the 1.1.1.1... in memory and start the next mix at 2.1.1.1.1... or just delete it?
        for example, let's say we split something into multiple streams in distillation column, there will be multiple streams with 1 2 3 4 5 etc, where recursive splits would be denoted as 2.1 2.2 etc.
        


    challenges:
        how do we prevent fluids from merging twice 
        (example: let's say that channel 1.1 merges with 2.2, when we go back to 2.2, how do we know that 2.2 has already been done?)
        (maybe just delete 2.2 and the solver will skip over it next time)
        
        how do we elegantly keep track of recursive splitting without it being a jumble of lists within lists 
        
        
        what would the ideal inputs look like?
            1. something like this? [[0, 1, 2, 3, 4], "split", [[5, 6, 7, 8, 9], [10, 11, 12, 13, 14], "mix"], [15, 16, 17, 18, 19]]
                no
            2. something like this? [[0, 1, 2, 3, 4], "split1", [5, 6, 7, 8, 9], "split1", [10, 11, 12, 13, 14], "mix1", [15, 16, 17, 18, 19]]
                maybe 
                it doesn't need a "mix1" after every "split1" stream, but it's impossible to tell which fluids should be mixed 
                but what happens if we have something like this, though: [[0, 1, 2, 3, 4], "split1", [5, 6, 7, 8, 9], "product", "split1", [10, 11, 12, 13, 14], "mix1"]
                how do we make sure that doesn't go back to the stream
            3. something like this? [[0, 1, 2, 3, 4], "split1", [5, 6, 7, 8, 9], "mix1", "split1", [10, 11, 12, 13, 14], "mix1", [15, 16, 17, 18, 19]]
                maybe (but better)
                the only issue would be not handling the "mix1" before the other "mix1"
                I guess for this one we could do this to show a product [[0, 1, 2, 3, 4], "split1", [5, 6, 7, 8, 9], "product", "split1", [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]]
            4. something like this? [[0, 1, 2, 3, 4], "split1" [[5, 6, 7, 8, 9] "mix1", [10, 11, 12, 13, 14] "mix1"], [15, 16, 17, 18, 19]]
                [[0, 1, 2, 3, 4], "split1" [[5, 6, 7, 8, 9] "product", [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]]]
            5. something like this? [[0, 1, 2, 3, 4], "split1" [[5, 6, 7, 8, 9], [10, 11, 12, 13, 14], "mix(0, 1)"], [15, 16, 17, 18, 19]]
                [[0, 1, 2, 3, 4], "split1" [[5, 6, 7, 8, 9], "product(0)", [10, 11, 12, 13, 14], "mix(1)"], [15, 16, 17, 18, 19]]
            6. something like this?
                [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14], [15, 16, 17, 18, 19]]
                [split(([0], [1, 2])), mix([1, 2], [3])]
              
                [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14], [15, 16, 17, 18, 19]]
                [split(([0], [1, 2])), product([1]))]
            7. something like this?
                [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14], [15, 16, 17, 18, 19]]
                [fluid1, fluid2, fluid3, fluid4]
                fluid2 and fluid3 are split mass fractions of fluid1 that are defined while this list is created
              
                [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14], [15, 16, 17, 18, 19]]
                [fluid1, fluid2, fluid3, fluid4]
            8. something like this?
                [[fluid1:[0, 1, 2, 3, 4], fluid2:[5, 6, 7, 8, 9], fluid3:[10, 11, 12, 13, 14], fluid4:[15, 16, 17, 18, 19]]
                fluid2 and fluid3 are split mass fractions of fluid1 that are defined while this list is created
                only issue is that we cannot define fluid4 (which would be a mix between fluid2 and fluid 3) before those fluids are run
              
                [[fluid1:[0, 1, 2, 3, 4], fluid2:[5, 6, 7, 8, 9], fluid3:[10, 11, 12, 13, 14], fluid4:[15, 16, 17, 18, 19]]
                only issue is that we cannot define fluid4 (which would be a mix between fluid2 and fluid3) before those fluids are run
                    this is made even worse when fluid2 or fluid3 are made product streams
                    this could be solved where each fluid is determined while the solver is running, although this would make too similar to 3, 4, 7, or 8
                    I still like this because it solves the issue of keeping track of fluid
                    another issue is that if there's a system where components in fluid2 system are run with a fluid for 1 hour and then a different one for 1 hour
                    however, I don't think most ChemE's are going to be using the same components for different fluids
            9. something like this?
                [[fluid0:[0, 1, 2, 3, 4], split(comp4, 1, 2), fluid1:[5, 6, 7, 8, 9], fluid2:[10, 11, 12, 13, 14], mix(1, 2), fluid3:[15, 16, 17, 18, 19]]
                or
                [[fluid0:[0, 1, 2, 3, 4], split(None, 1, 2), fluid1:[5, 6, 7, 8, 9], fluid2:[10, 11, 12, 13, 14], mix(1, 2), fluid3:[15, 16, 17, 18, 19]]
                where the split checks if component 4 is a splitter, and if it isn't it just uses a 0.5 / 0.5 mass split

                [[fluid0:[0, 1, 2, 3, 4], split(None, 1, 2), fluid1:[5, 6, 7, 8, 9], product(1), fluid2:[10, 11, 12, 13, 14, 15, 16, 17, 18, 19]]
                for products
                
                for heat exchangers and other things:
                    heat_exchange(comp1, comp10)
                [[fluid0:[0, 1, 2, 3, 4], split(None, 1, 2), fluid1:[5, 6, 7, 8, 9], fluid2:[10, 11, 12, 13, 14], heat_exchange(comp2, comp10), mix(1, 2), fluid3:[15, 16, 17, 18, 19]]
            10. something like this?
                this one uses component values instead of dictionary values for mixers and splitters
                [[fluid0:[0, 1, 2, 3, 4], split(4, 5, 10), fluid1:[5, 6, 7, 8, 9], fluid2:[10, 11, 12, 13, 14], mix(9, 14), fluid3:[15, 16, 17, 18, 19]]
                where the split checks if component 4 is a splitter, and if it isn't it just uses a 0.5 / 0.5 mass split

                [[fluid0:[0, 1, 2, 3, 4], split(4, 5, 10), fluid1:[5, 6, 7, 8, 9], product(9), fluid2:[10, 11, 12, 13, 14, 15, 16, 17, 18, 19]]
                for products
                
                for heat exchangers and other things:
                heat_exchange(1, 10)
                [[fluid0:[0, 1, 2, 3, 4], split(4, 5, 10), fluid1:[5, 6, 7, 8, 9], product(9), fluid2:[10, 11, 12, 13, 14, 15, 16, 17, 18, 19]]
                we could either do this for heat exchange or whenever a heat exchange object is encountered, the program uses that object to do the calculations
            
            11. something like this?
                this one will just have a basic array and there will be if statements in the loop that will determine how to handle splits and mixing
                
            I think my favourite is either 7, 9, or 10
"""
        
        
    
        
            
    #def clear_results():

#wow this is really clean
zero_compressor = compressor(None, "compressor 0", None, 3, 0.005, 0.007, 0.65, 100000, 10000, 0, 0)

first_compressor = compressor([zero_compressor], "compressor 1", None, 3, 0.005, 0.007, 0.65, -1000, 0, 0, 0)
second_compressor = compressor([zero_compressor], "compressor 2", None, 3, 0.005, 0.007, 0.65, 0, 100000, 0, 0)

third_compressor = compressor([first_compressor, second_compressor], "compressor 3", None, 3, 0.005, 0.007, 0.65, 50000, 0, 0, 0)
fourth_compressor = compressor([third_compressor], "compressor 4", None, 3, 0.005, 0.007, 0.65, 50000, 0, 0, 0)
zero_compressor.upstreams = [fourth_compressor]


print("test", third_compressor.work_in)

print(first_compressor.enthalpy_change)


#FIX ME: using downstreams here is problematic, while it's a lot easier to keep track of, variables cannot be accessed before they are initialized :/

#before_compressor_test2.downstreams = [compressor2]


starting_stream = ct.Solution("gri30.yaml", "gri30")
starting_stream.TPX = 400, 100000, "H2O:1, CH4:1"
#print(starting_stream.h)


Rankine_Cycle = system(starting_stream, zero_compressor, [zero_compressor, first_compressor, second_compressor, third_compressor, fourth_compressor])
#Rankine_Cycle.process_comp_list()
#Rankine_Cycle.run()
Rankine_Cycle.better_run()


#manager.show_streams()
#print(starting_stream.TPX)
#print(starting_stream.h)
#print(compressor2.enthalpy_change)
#fluid_check = manager.get_stream(current_stream)

def mix_streams_with_upstreams(input_component):
    """
        example input:
        split_streams(splitter1, {compressor3: 0.5, compressor9: 0.5})
        split_streams(pipe1, {pipe2: 0.25, pipe3: 0.75})
        
        here's the referenced cantera docs:
        https://cantera.org/3.1/examples/python/reactors/mix1.html 
    """
    reservoirs = []
    
    total_mass_flow = sum(input_component.upstreams[i].mass_flow for i in range(len(input_component.upstreams)))

    for i in range(len(input_component.upstreams)):
        reservoirs.append(ct.Reservoir(input_component.upstreams[i].outlet_fluids, name=None))
    
    mixer = ct.IdealGasReactor(input_component.inlet_fluids[0], name=None)
    
    for i in range(len(input_component.inlet_fluids)):
        ct.MassFlowController(reservoirs[i], mixer, input_component.upstreams[i].mass_flow)
        
    output_reservoir = ct.Reservoir(input_component.inlet_fluids[0])
    
    #high K value to not hinder flow 
    outlet = ct.Valve(mixer, output_reservoir, K=1e9)
    
    sim = ct.ReactorNet([mixer])
    
    sim.advance_to_steady_state()
    
    input_component.outlet_fluids = mixer.thermo
    input_component.mass_flow = total_mass_flow
    


