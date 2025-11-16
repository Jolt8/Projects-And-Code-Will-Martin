import numpy as np

import math 

from scipy import optimize

import ht as ht

#from HeatExchangers_Thermo_Final import heat_transfer_coeff

import networkx as nx

from thermo import ChemicalConstantsPackage, CEOSGas, PRMIX, CEOSLiquid, FlashVL, equilibrium, EquilibriumStream, EquilibriumState, eos_mix, Flash, Stream, StreamArgs

from enum import Enum, auto



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
        
        
        Compressor 
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
"""


""" Energy Related Classes
- An ideal equalizer (ideal heat exchanger)
- A non-ideal equalizer (non ideal heat exchanger)


For adding/removing energy, we could either use 
- A energy creator class and destroyer that we link to a forced energy stream class
- Each energy stream class has options to input either a negative or positive wattage value 

For transferring heat with fluids we could either use:
- Dedicated heat moving streams with fluids 
- Each energy stream class has the option to specify whether it uses fluids to transfer 
    - Think like heat pipes, circulated water for cooling, etc.
- This would return the required fluid flow for moving a set amount of heat or equalizing two components to within a certain toleranrce
- NOTE: This tolerance idea for equalizing is essential, don't forget

Should we have different classes for different energy types, or should we just specify what type of energy they are?
- for heat, electrical, mechanical, etc. 

We should probably specify these types of parameters (just writing down ideas)
- transmission efficiency - either set or calculated due to heat loss
    - Maybe wire guage, voltage, amperage, frequency could be specified to calculate this for electrical systems
    - Not sure what we could do to calculate friction losses for mechanical systems, we'll probably just use either a percentage loss or a set enthalpy loss value
    - in each of these cases, is there a good way to handle three sets of input methods that all take different amounts of inputs?
        - For example, for the electrical one outlined above, having a simple percent losses, another static enthalpy losses, and a third 
        method for the complicated wire guage, voltage, amperage, etc.
        
        
How do we handle components where they have inherent wattage inputs 
- For example, for a distillation column object that has multiple wattage inputs for each stage (mostly for reboiler and condenser though), 
should we have the distillation column require inputs from energy streams, have the user specify the column's wattage, or a hybrid approach

Currently, the way I define my components is through specifying their upstreams, should I just specify both the upstreams and downstreams for energy streams

Should the way I handle energy streams be pretty similar to how I handle components in terms of looping logic 
"""

"""Some Additional Questions:

Should there be two seperate loops where one is for components and then the other is for energy streams?

Should I use a big networkx graph that contains both components and energy streams?

What do we do in instances where the wattage of an upstream component is sent to a downstream one and the cycle (for components and not energys streams at least) 
is not cyclical (We'll probably just identify this all the same and wait until the system converges into a steady state)
"""




"""Future Things we Should Probably Add:

A way to find components online that match or closely match your specificiations 
    - inputting the link to a component and using its specs within the simulator 
    - something like KG_Tower but more universal
    - Maybe I should just focus on this instead of the process sim and start developing this for aspen plus or HYSYS
"""

class EnergyType(Enum):
    thermal = auto()
    mechanical = auto()
    electrical = auto()

class LossModel():
    def __init__(self):
        pass


class Percentage():
    def __init__(self, percentage_loss):
        LossModel.__init__(self)
   

class Fixed():
    def __init__(self, fixed_wattage_loss):
        LossModel.__init__(self)
        

class DetailedElectrical():
    def __init__(self, wire_guage, length, voltage, amperage, frequency):
        LossModel.__init__(self)
        

class Frictional():
    def __init__(self, friction_coeff:list):
        LossModel.__init__(self)

class Component_Template():
    def __init__(self, upstreams, name):
        self.upstreams = upstreams if upstreams is None else upstreams
        self.name = name
        self.inlet_streams = []
        
        
class EnergyStream():
    def __init__(self, power, energy_type:EnergyType, ):
        self.power = power
        self.energy_type = energy_type

class EnergySource(Component_Template):
    def __init__(self, name, power, energy_type:EnergyType, loss_model:LossModel):
        super().__init__(None, name)
        self.outlet_streams = []
        self.power = power
        self.energy_type = energy_type
        self.loss_model = loss_model
        
    def update_state(self):
        self.outlet_streams.append(EnergyStream(self.power, self.energy_type))
        
    
        if self.loss_model == Percentage:
            self.power = self.power * (1 - self.loss_model.percentage_loss)
            
        elif self.loss_model == Fixed:
            self.power = self.power - self.loss_model.fixed_loss
            
        elif self.loss_model == DetailedElectrical:
            #example
            power_loss = (self.power / self.loss_model.wire_guage * self.loss_model.length)
            self.power = self.power - power_loss
            
        elif self.loss_model == Frictional:
            #example
            power_loss = (self.power / sum(self.loss_model.friction_coeff))
            self.power = self.power - power_loss
            
        else:
            pass
    

class Pipe(Component_Template):
    def __init__(self, upstreams:list, name, inlet_streams:list, molar_flow, inner_characteristic_length, outer_characteristic_length, 
                 length, do_heat_loss=False, pipe_thermal_conductivity=0, expected_temperature_rise=30, air_temperature_initial=300, air_temperature_final=320, 
                 inner_wall_temperature=320, outer_wall_temperature=280, 
                 do_pressure_drop=False, pressure_drop_method="simple", absolute_pipe_roughness=0, 
                 ):
        """
        For pressure drop method, your two options are "simple" and "twophase"
        """
        Component_Template.__init__(self, upstreams, name, inlet_streams, molar_flow)
        self.outlet_streams = []
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
                self.inlet_streams[0], self.molar_flow, self.inner_characteristic_length, self.length, self.absolute_pipe_roughness)
            self.current_fluid.P = self.current_fluid.P - self.pressure_drop
        
        elif self.pressure_drop_method == "twophase":
            self.pressure_drop = cantera_incremental_pressure_drop(self.inlet_streams[0], self.molar_flow, self.inner_characteristic_length, self.length, self.absolute_pipe_roughness)
            self.current_fluid.P = self.current_fluid.P - self.pressure_drop
            
    def heat_change(self):
        #FIXME: I need a good way to find the final temperatures of the internal fluid and the rise in air temperature
        self.heat_exchange = cantera_air_heat_loss(self.inlet_streams[0], self.expected_temperature_rise, self.molar_flow, 
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
            
            
class Splitter(Component_Template):
    """
    Example of outlet_stream_split variable: [0.8, 0.8, 9]
    """
    def __init__(self, upstreams:list, name, inlet_streams:list, inner_characteristic_length, outer_characteristic_length, respective_split_ratios:list):
        Component_Template.__init__(self, upstreams, name, inlet_streams, inner_characteristic_length, outer_characteristic_length)
        self.outlet_streams = []
        
        
class Connect():
    def __init__(self, comp_from, comp_to):
        if comp_from not in comp_to.upstreams:
            comp_to.upstreams.append(comp_from)
        else:
            print(comp_from, "is already connected to", comp_to)


class Valve(Component_Template):
    #Should this be for both valves that drop pressure and limit flow 
    def __init__(self, upstreams:list, name, inner_characteristic_length, outer_characteristic_length, pressure_drop, pressure_drop_ratio, molar_flow_limit, ):
        Component_Template.__init__(self, upstreams, name, inner_characteristic_length, outer_characteristic_length)
        self.outlet_streams = []
        self.pressure_drop = pressure_drop
        self.pressure_drop_ratio = pressure_drop_ratio
        self.molar_flow_limit = molar_flow_limit
    
    def update_state(self):
        pass
        #oooh noooo, whenever this limits flow we have to look ahead in the simulator and adjust the flow rate for other components 


class GasTurbine(Component_Template):
    def __init__(self, upstreams:list, name, inner_characteristic_length, outer_characteristic_length, efficiency, pressure_change, enthalpy_change, work_in, compression_ratio):
        Component_Template.__init__(self, upstreams, name, inner_characteristic_length, outer_characteristic_length)
        self.outlet_streams = []
        self.efficiency = efficiency


class Compressor(Component_Template):
    def __init__(self, upstreams:list, name, efficiency, pressure_change, enthalpy_change, compression_ratio):
        Component_Template.__init__(self, upstreams, name)
        self.inlet_streams = []
        self.inlet_energy_streams = []
        self.outlet_streams = []
        self.mat_upstreams = []
        self.e_upstreams = []
        self.efficiency = efficiency
        self.pressure_change = pressure_change
        self.enthalpy_change = enthalpy_change
        self.compression_ratio = compression_ratio
        self.permissable_energy_streams = [EnergyType.electrical, EnergyType.mechanical]
    
    
    def update_state(self):
        outlet_stream = Stream(IDs=self.inlet_streams[0].IDs, T=self.inlet_streams[0].T, P=self.inlet_streams[0].P, zs=self.inlet_streams[0].zs, n=self.inlet_streams[0].n)
        
        if self.inlet_energy_streams != []:
            if self.inlet_energy_streams[0].energy_type not in self.permissable_energy_streams:
                (print("Wrong energy source", self.inlet_energy_streams.energy_type))
        
        if self.pressure_change != 0:
            final_pressure = self.inlet_streams[0].P + self.pressure_change
            outlet_stream.P = final_pressure
            
            
            state_2_ideal = Stream(IDs=self.inlet_streams[0].IDs, S=self.inlet_streams[0].S, P=final_pressure, zs=self.inlet_streams[0].zs, n=self.inlet_streams[0].n)
            delta_H_ideal = (state_2_ideal.H - self.inlet_streams[0].H)
            # The definition of isentropic efficiency means that the actual amount of heat added is
            # dH_actual = dH_ideal/eta_isentropic
            H_added_to_fluid_actual = delta_H_ideal / self.efficiency
            outlet_stream.flash(H=self.inlet_streams[0].H + H_added_to_fluid_actual, P=final_pressure)
        #the enthalpy change one might be unnecessary 
        elif self.enthalpy_change != 0: #^ self.enthalpy_change is not None:
            #TODO: his assumes constant pressure, make sure to fix this later
            outlet_stream.flash(H=(self.inlet_streams[0].H + self.enthalpy_change), P=self.inlet_streams[0].P)
        elif self.inlet_energy_streams != []: 
            initial_pressure = outlet_stream.P
            initial_enthalpy = outlet_stream.H
            
            self.enthalpy_change = self.inlet_energy_streams[0].power / self.inlet_streams[0].n
            
            new_enthalpy = initial_enthalpy + (self.enthalpy_change / self.efficiency)
            outlet_stream.H = new_enthalpy
            
            self.pressure_change = initial_pressure - outlet_stream.P
        elif self.compression_ratio != 0: #^ self.compression_ratio is not None:
            final_pressure = self.inlet_streams[0].P * self.compression_ratio
            state_2_ideal = Stream(IDs=self.inlet_streams[0].IDs, S=self.inlet_streams[0].S, P=final_pressure, zs=self.inlet_streams[0].zs, n=self.inlet_streams[0].n)
            
            # Compute the change in enthalpy
            delta_H_ideal = (state_2_ideal.H-self.inlet_streams[0].H)
            # The definition of isentropic efficiency means that the actual amount of heat added is
            # dH_actual = dH_ideal/eta_isentropic
            H_added_to_fluid_actual = delta_H_ideal * self.efficiency
            outlet_stream.flash(H=self.inlet_streams[0].H + H_added_to_fluid_actual, P=final_pressure)
        self.outlet_streams.append(outlet_stream)

def mix_mat_streams_with_upstreams(input_component):
    
    total_moles = 0
    for i in range(len(input_component.inlet_streams)):
        total_moles += input_component.inlet_streams[i].n
    
    """
    for i in range(len(input_component.inlet_streams)):
        print("hey pal", input_component.inlet_streams[i].H)
    """
    
    
    total_enthalpy = sum(fluid.H * fluid.n for fluid in input_component.inlet_streams)

    zs_mixed = [sum(fluid.zs[i] * fluid.n for fluid in input_component.inlet_streams) / total_moles for i in range(len(input_component.inlet_streams[0].zs))]

    # You need to decide on the outlet pressure. Often, it's the lowest inlet pressure.
    mixed_pressure = min(fluid.P for fluid in input_component.inlet_streams)

    # Perform an enthalpy-pressure (HP) flash
    mixed_state = Stream(IDs=input_component.inlet_streams[0].IDs, zs=zs_mixed, H=(total_enthalpy/total_moles), P=mixed_pressure, n=input_component.inlet_streams[i].n) # Note H is molar here
    input_component.inlet_streams = [mixed_state]


def mix_e_streams_with_upstreams(input_component):
    total_power = 0
    for i in range(len(input_component.inlet_energy_streams)):
        total_power += input_component.inlet_energy_streams[i].power
        
    input_component.inlet_streams = EnergyStream(total_power, EnergyType.thermal, None)

class System():
    def __init__(self, start_state, starting_component, components):
        self.start_state = start_state
        self.starting_component = starting_component
        self.components = components
    
    def get_total_head_loss(self):
        total_head_loss = []
        for i in range(len(self.components)):
            if hasattr(self.components[i], "Pipe"):
                total_head_loss.append(self.components[i].head_loss)
                 
    
    def run(self):
        self.get_total_head_loss = System.get_total_head_loss(self)
        #order the list starting from the starting position defined by the user
        self.components = self.components[self.components.index(self.starting_component):] + self.components[:self.components.index(self.starting_component)]
        
        G = nx.DiGraph()
        G.add_nodes_from(self.components)
        
        #you have to do this or else you will get an error that the dictionary changed during iteration 
        nodes = list(G.nodes)
        for node in nodes:
            connections = []
            #this should probably be be edited later so that it checks if something is a dead end or start or not 
            if isinstance(node, EnergySource):
                pass
            else:
                for i in range(len(node.upstreams)):
                    connections.append([node.upstreams[i], node])
            G.add_edges_from((connections))
        
            #G.add_edges_from([(node.upstreams[ds_index], node) for ds_index in range(len(node.upstreams))])
       
        visited = [] #change back to set() when done
        
        
        #print([self.components[i].name for i in range(len(self.components))])
        #print("hey pal", self.starting_component.inlet_streams, self.start_state)
        self.starting_component.inlet_streams = [self.start_state]
        #print("hey bud", self.starting_component.inlet_streams, self.start_state)
        self.starting_component.update_state() 
        visited.append(self.starting_component.name)
        
        for i in range(len(self.components)):
            mat_upstreams = []
            e_upstreams = []
            if self.components[i].upstreams != None:
                for j in range(len(self.components[i].upstreams)):
                    if hasattr(self.components[i].upstreams[j], "energy_type"):
                        print("success")
                        e_upstreams.append(self.components[i].upstreams[j])
                    else:
                        mat_upstreams.append(self.components[i].upstreams[j])
                self.components[i].mat_upstreams = mat_upstreams
                self.components[i].e_upstreams = e_upstreams
       
        print(f"\n--- {self.starting_component.name} ---")
        print(f"Inlet: T={self.starting_component.inlet_streams[0].T:.2f}K, P={self.starting_component.inlet_streams[0].P/1000:.2f}kPa, n={self.starting_component.inlet_streams[0].n:.2f}mol, H={self.starting_component.inlet_streams[0].H/1000:.2f}kJ/mol")
        print(f"Outlet: T={self.starting_component.outlet_streams[0].T:.2f}K, P={self.starting_component.outlet_streams[0].P/1000:.2f}kPa, n={self.starting_component.outlet_streams[0].n:.2f}mol, H={self.starting_component.outlet_streams[0].H/1000:.2f}kJ/mol")
        
        start_index = self.components.index(self.starting_component)
        
        for src, node in nx.edge_bfs(G, self.components[start_index]):
            if node.name not in visited: 
                
                print(f"\n--- {node.name} ---")
                #This handles basically everything including splitting, 
                total_moles = 0
                for i in range(len(node.mat_upstreams)):
                    for j in range(len(node.mat_upstreams[i].outlet_streams)):
                        total_moles += node.mat_upstreams[i].outlet_streams[j].n / G.out_degree(node.mat_upstreams[i])
                for i in range(len(node.mat_upstreams)):
                    print(node.mat_upstreams[i].outlet_streams)
                    for j in range(len(node.mat_upstreams[i].outlet_streams)):
                        node.inlet_streams.append(Stream(IDs=node.mat_upstreams[i].outlet_streams[j].IDs, T=node.mat_upstreams[i].outlet_streams[j].T, 
                                                        P=node.mat_upstreams[i].outlet_streams[j].P, zs=node.mat_upstreams[i].outlet_streams[j].zs, 
                                                        n=total_moles))
                total_power = 0
                for i in range(len(node.e_upstreams)):
                    for j in range(len(node.e_upstreams[i].outlet_streams)):
                        total_power += node.e_upstreams[i].outlet_streams[j].power / G.out_degree(node.e_upstreams[i])
                for i in range(len(node.e_upstreams)):
                    print(node.e_upstreams[i].outlet_streams)
                    if node.e_upstreams[i].upstreams == None:
                            node.e_upstreams[i].update_state()
                    for j in range(len(node.e_upstreams[i].outlet_streams)):
                        node.inlet_energy_streams.append(node.e_upstreams[i].outlet_streams[j])
                    
                
                #this is for mixing states (the current component has two inputs)
                if len(node.mat_upstreams) >= 2:
                    mix_mat_streams_with_upstreams(node)
                    
                if len(node.e_upstreams) >= 2:
                    mix_e_streams_with_upstreams(node)
                
                print(node.name)
                node.update_state()
                visited.append(node.name)
                
                print(f"Inlet: T={node.inlet_streams[0].T:.2f}K, P={node.inlet_streams[0].P/1000:.2f}kPa, n={node.inlet_streams[0].n:.2f}mol, H={node.inlet_streams[0].H/1000:.2f}kJ/mol")
                print(f"Outlet: T={node.outlet_streams[0].T:.2f}K, P={node.outlet_streams[0].P/1000:.2f}kPa, n={node.outlet_streams[0].n:.2f}mol, H={node.outlet_streams[0].H/1000:.2f}kJ/mol")
                if node.name == "compressor 4":
                    print(node.inlet_energy_streams[0].power)

                
        new_dict = []
        for src, node in nx.edge_bfs(G, self.components[start_index]):
            new_dict.append((src.name, "-->", node.name))
        print(new_dict)

        
            
    #def clear_results():

#wow this is really clean
#mass flow for only the input component or if you want to add more fluid (although you probably want to define the state first)
zero_compressor = Compressor(None, "compressor 0", 0.65, 10000, 0, 0)

first_compressor = Compressor([zero_compressor], "compressor 1", 0.65, 0, 0, 10)
second_compressor = Compressor([zero_compressor], "compressor 2", 0.65, 0, 1000000, 0)

third_compressor = Compressor([first_compressor, second_compressor], "compressor 3", 0.65, 50000, 0, 0)

zero_energy_source = EnergySource("Energy Source 0", 5000, EnergyType.electrical, None)
fourth_compressor = Compressor([third_compressor, zero_energy_source], "compressor 4", 0.65, 0, 0, 0)

zero_compressor.upstreams = [fourth_compressor]

zero_wall = 1

#print(first_compressor.enthalpy_change)

#FIXME: using dowstreams here is problematic, while it's a lot easier to keep track of, variables cannot be accessed before they are initialized :/
#This was fixed, however I still feel like the thinking in terms of debugging, working with these is sometimes problematic, 
#however remember that this doesn't apply to how the user defines stuff, as shown above

#before_compressor_test2.downstreams = [compressor2]

IDs = ['Water', 'Carbon Dioxide']
constants, correlations = ChemicalConstantsPackage.from_IDs(IDs=IDs)
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)
stream = Stream(IDs=IDs, T=350, P=100000, zs=[0.5, 0.5], n=3)


Rankine_Cycle = System(stream, zero_compressor, [zero_compressor, first_compressor, second_compressor, third_compressor, zero_energy_source, fourth_compressor])
Rankine_Cycle.run()

#print(zero_compressor.outlet_streams)


#this is just for future convergence logic
"""
for i in range(5):
    new_stream = Stream(IDs=fourth_compressor.outlet_streams[0].IDs, T=fourth_compressor.outlet_streams[0].T, 
                        P=fourth_compressor.outlet_streams[0].P, zs=fourth_compressor.outlet_streams[0].zs, 
                        n=fourth_compressor.outlet_streams[0].n)
    
    
    
    print("test", zero_compressor.inlet_streams)
    new_cycle = System(new_stream, zero_compressor, [zero_compressor, first_compressor, second_compressor, third_compressor, fourth_compressor])
    for i in range(len(new_cycle.components)):
        new_cycle.components[i].outlet_streams = []
        new_cycle.components[i].inlet_streams = []
    zero_compressor.inlet_streams = [new_stream]
    new_cycle.run()
"""

#print(second_compressor.inlet_streams[0].T, second_compressor.outlet_streams)