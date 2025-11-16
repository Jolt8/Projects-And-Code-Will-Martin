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

#from HeatExchangers_Thermo_Final import heat_transfer_coeff

import networkx as nx

from thermo import ChemicalConstantsPackage, CEOSGas, PRMIX, CEOSLiquid, FlashVL, equilibrium, EquilibriumStream, EquilibriumState, eos_mix, Flash



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
        


"""
print("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEE", input_component.inlet_fluids)
print(input_component.inlet_fluids[0].TPX)
print(input_component.inlet_fluids[1].TPX)


A = ct.Quantity(input_component.inlet_fluids[0], constant="HP")
quantities = ct.Quantity((input_component.input_components[i] for i in range(len(input_component.inlet_fluids))), constant="HP", )
input_component.inlet_fluids = input_component.inlet_fluids[0] + input_component.inlet_fluids[1]
input_component.molar_flow = mixer.molar_flow_rates
"""


    

    



class component_template():
    def __init__(self, upstreams:list, name, inlet_fluids:list, molar_flow, inner_characteristic_length, outer_characteristic_length):
        #how do I make it so that this component doesn't require stream inputs from the user but 
        self.upstreams = [] if upstreams is None else upstreams
        self.name = name
        self.inlet_fluids = [] if inlet_fluids is None else inlet_fluids
        self.molar_flow = molar_flow
        #self.split_count = 0
        self.inner_characteristic_length = inner_characteristic_length
        self.outer_characteristic_length = outer_characteristic_length
        self.thickness = self.inner_characteristic_length - self.outer_characteristic_length
        
        try:
            for i in range(len(self.inlet_fluids)):
                density = inlet_fluids[i].density
                
                area = (math.pi * ((inner_characteristic_length / 2) ** 2))
                
                if self.molar_flow is not None:
                    self.volumetric_flow = molar_flow / density 
                    self.velocity = self.volumetric_flow / area
                elif self.volumetric_flow is not None:
                    self.molar_flow = self.volumetric_flow / density
                    self.velocity = self.volumetric_flow / area
                elif self.velocity is not None:
                    self.volumetric_flow = area / self.velocity
                    self.molar_flow = self.volumetric_flow / density
        except:
            pass
            
        #def environmental_heat_exchange():
            


class Expander(component_template):
    def __init__(self, upstreams:list, name, inlet_fluids:list, molar_flow, inner_characteristic_length, outer_characteristic_length, efficiency, pressure_change, enthalpy_change, work_in, compression_ratio):
        component_template.__init__(self, upstreams, name, inlet_fluids, molar_flow, inner_characteristic_length, outer_characteristic_length)
        self.efficiency = efficiency
        self.pressure_change = pressure_change
        self.enthalpy_change = enthalpy_change
        self.work_in = work_in
        self.compression_ratio = compression_ratio
        self.outlet_fluids = []
    
    def update_state(self):
        
        outlet_fluid = ct.Solution("gri30.yaml", "gri30")
        outlet_fluid.TPX = self.inlet_fluids[0].TPX
        
        if self.pressure_change != 0:
            initial_temperature = outlet_fluid.T
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            outlet_fluid.SP = initial_entropy, (initial_pressure + self.pressure_change)
            new_enthalpy = outlet_fluid.h
            
            isentropic_enthalpy_change = new_enthalpy - initial_enthalpy
            self.enthalpy_change = isentropic_enthalpy_change / self.efficiency
            self.work_in = self.enthalpy_change * self.molar_flow
            outlet_fluid.HP = initial_enthalpy + self.enthalpy_change, (initial_pressure + self.pressure_change)
        #the enthalpy change one might be unnecessary 
        elif self.enthalpy_change != 0: #self.enthalpy_change is not None:
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            new_enthalpy = initial_enthalpy + self.enthalpy_change
            outlet_fluid.HP = new_enthalpy, None
            
            self.pressure_change = initial_pressure - outlet_fluid.P
            self.work_in = self.enthalpy_change * self.molar_flow
        elif self.work_in != 0: #self.work_in is not None:
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            self.enthalpy_change = self.work_in / self.molar_flow
            
            new_enthalpy = initial_enthalpy + self.enthalpy_change
            outlet_fluid.HP = new_enthalpy, None
            
            self.pressure_change = initial_pressure - outlet_fluid.P
        elif self.compression_ratio != 0: #self.compression_ratio is not None:
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.s
            
            self.pressure_change = initial_pressure * self.compression_ratio - initial_pressure
            
            outlet_fluid.SP = initial_entropy, (initial_pressure + self.pressure_change)
            new_enthalpy = outlet_fluid.h
            
            isentropic_enthalpy_change = new_enthalpy - initial_enthalpy
            self.enthalpy_change = isentropic_enthalpy_change / self.efficiency
            self.work_in = self.enthalpy_change * self.molar_flow
            outlet_fluid.HP = initial_enthalpy + self.enthalpy_change, (initial_pressure + self.pressure_change)
        self.outlet_fluids.append(outlet_fluid)

class Pipe(component_template):
    def __init__(self, upstreams:list, name, inlet_fluids:list, molar_flow, inner_characteristic_length, outer_characteristic_length, 
                 length, do_heat_loss=False, pipe_thermal_conductivity=0, expected_temperature_rise=30, air_temperature_initial=300, air_temperature_final=320, 
                 inner_wall_temperature=320, outer_wall_temperature=280, 
                 do_pressure_drop=False, pressure_drop_method="simple", absolute_pipe_roughness=0, 
                 ):
        """For pressure drop method, your two options are "simple" and "twophase"
        
        
        """
        component_template.__init__(self, upstreams, name, inlet_fluids, molar_flow, inner_characteristic_length, outer_characteristic_length)
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
                self.inlet_fluids[0], self.molar_flow, self.inner_characteristic_length, self.length, self.absolute_pipe_roughness)
            self.current_fluid.P = self.current_fluid.P - self.pressure_drop
        
        elif self.pressure_drop_method == "twophase":
            self.pressure_drop = cantera_incremental_pressure_drop(self.inlet_fluids[0], self.molar_flow, self.inner_characteristic_length, self.length, self.absolute_pipe_roughness)
            self.current_fluid.P = self.current_fluid.P - self.pressure_drop
            
    def heat_change(self):
        #FIXME: I need a good way to find the final temperatures of the internal fluid and the rise in air temperature
        self.heat_exchange = cantera_air_heat_loss(self.inlet_fluids[0], self.expected_temperature_rise, self.molar_flow, 
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
    def __init__(self, upstreams:list, name, inlet_fluids:list, molar_flow, inner_characteristic_length, outer_characteristic_length, respective_split_ratios:list):
        component_template.__init__(self, upstreams, name, inlet_fluids, molar_flow, inner_characteristic_length, outer_characteristic_length)
        self.outlet_fluids = []
        
        

"""
    oh shit! how do we do heat / energy transfer?
        walls
            wall1(adiabadic_reactor, pipe1, wall_effectiveness)
                wall_effectiveness describes the ratio between the heat transfered into the other component and heat still lost to environment
                    ex. if a Compressor or reactor was creating 500w of heat, a ratio of 0.5 would allow 250w to be actually transfered, the other 250w lost
                    if we didn't do this at all, we would lose 500w of heat regardless
        energy_transfer_fluid
            energy_transfer_fluid(adiabadic_reactor, pipe1, )
                for heat transfer using a working fluid
        wire / axle
            wire / axle(expander, pump)
"""

class Compressor(component_template):
    def __init__(self, upstreams:list, name, inlet_fluids:list, molar_flow, inner_characteristic_length, outer_characteristic_length, efficiency, pressure_change, enthalpy_change, work_in, compression_ratio):
        component_template.__init__(self, upstreams, name, inlet_fluids, molar_flow, inner_characteristic_length, outer_characteristic_length)
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
        
        outlet_fluid = flasher.flash(self.inlet_fluids[0].zs, self.inlet_fluids[0].T, self.inlet_fluids[0].P)
        
        if self.pressure_change != 0:
            print("pressure_change")
            final_pressure = self.inlet_fluids[0].P + self.pressure_change
            state_2_ideal = flasher.flash(S=self.inlet_fluids[0].S(), P=final_pressure, zs=self.inlet_fluids[0].zs)
            # Compute the change in enthalpy
            delta_H_ideal = (state_2_ideal.H()-self.inlet_fluids[0].H())
            # The definition of isentropic efficiency means that the actual amount of heat added is
            # dH_actual = dH_ideal/eta_isentropic
            H_added_to_fluid_actual = delta_H_ideal/self.efficiency
            outlet_fluid = flasher.flash(H=self.inlet_fluids[0].H() + H_added_to_fluid_actual, P=final_pressure, zs=self.inlet_fluids[0].zs)
        #the enthalpy change one might be unnecessary 
        elif self.enthalpy_change != 0: #^ self.enthalpy_change is not None:
            print("enthalpy_change")
            #TODO: his assumes constant pressure, make sure to fix this later
            
            outlet_fluid = flasher.flash(U=(self.inlet_fluids[0].U() + self.enthalpy_change), P=self.inlet_fluids[0].P, zs=self.inlet_fluids[0].zs)
        elif self.work_in != 0: #^ self.work_in is not None:
            print("work_in")
            initial_pressure = outlet_fluid.P
            initial_enthalpy = outlet_fluid.h
            initial_entropy = outlet_fluid.S()
            
            self.enthalpy_change = self.work_in / self.molar_flow
            
            new_enthalpy = initial_enthalpy + self.enthalpy_change
            outlet_fluid.H = new_enthalpy
            
            self.pressure_change = initial_pressure - outlet_fluid.P
        elif self.compression_ratio != 0: #^ self.compression_ratio is not None:
            print("compression_ratio")
            final_pressure = self.inlet_fluids[0].P * self.compression_ratio
            state_2_ideal = flasher.flash(S=self.inlet_fluids[0].S(), P=final_pressure, zs=self.inlet_fluids[0].zs)
            # Compute the change in enthalpy
            delta_H_ideal = (state_2_ideal.H()-self.inlet_fluids[0].H())
            # The definition of isentropic efficiency means that the actual amount of heat added is
            # dH_actual = dH_ideal/eta_isentropic
            H_added_to_fluid_actual = delta_H_ideal/self.efficiency
            outlet_fluid = flasher.flash(H=self.inlet_fluids[0].H() + H_added_to_fluid_actual, P=final_pressure, zs=self.inlet_fluids[0].zs)
        self.outlet_fluids.append(outlet_fluid)

def mix_streams_with_upstreams(input_component):
    """
        example input:
        split_streams(splitter1, {compressor3: 0.5, compressor9: 0.5})
        split_streams(pipe1, {pipe2: 0.25, pipe3: 0.75})
    """
    """
    temps = []
    pressures = []
    zss = []
    
        temps.append(input_component.inlet_fluids[i].T)
        pressures.append(input_component.inlet_fluids[i].P)
        
    
    zipped_zss = [sum(items) / len(input_component.upstreams) for items in zip(*zss)]
        
    
    mix_temp = sum(temps) / len(input_component.upstreams)
    mix_pressure = sum(pressures) / len(input_component.upstreams)
        
    mixed_state = flasher.flash(zs=zipped_zss, T=mix_temp, P=mix_pressure)
    
    input_component.inlet_fluids[0] = mixed_state
    """
    
    total_moles = input_component.molar_flow
    total_enthalpy = sum(fluid.H() * input_component.molar_flow for fluid in input_component.inlet_fluids)
    
    zs_mixed = [sum(fluid.zs[i] * input_component.molar_flow for fluid in input_component.inlet_fluids) / total_moles for i in range(len(input_component.inlet_fluids[0].zs))]
    
    zs_normalized = []
    for i in range(len(zs_mixed)):
        zs_normalized.append(zs_mixed[i] / len(input_component.upstreams))

    print(total_moles, total_enthalpy, zs_normalized)
    # You need to decide on the outlet pressure. Often, it's the lowest inlet pressure.
    mixed_pressure = min(fluid.P for fluid in input_component.inlet_fluids)

    # Perform an enthalpy-pressure (HP) flash
    mixed_state = flasher.flash(zs=zs_normalized, H=total_enthalpy/total_moles, P=mixed_pressure) # Note H is molar here
    input_component.inlet_fluids = [mixed_state]

class System():
    def __init__(self, start_state, flasher, starting_component, components):
        self.start_state = start_state
        self.starting_component = starting_component
        self.components = components
        self.total_molar_flow = self.starting_component.molar_flow
        #we should probably implement something like this
        #self.total_heat_wattage = total_heat_wattage
    
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
            G.add_edges_from([(node.upstreams[ds_index], node) for ds_index in range(len(node.upstreams))])
       
        visited = [] #change back to set() when done
        
        #print([self.components[i].name for i in range(len(self.components))])
        self.starting_component.inlet_fluids = [self.start_state]
        self.starting_component.update_state() 
        visited.append(self.starting_component.name)
        
        print("")
        print(self.starting_component.name)
        print("current component previous connections:", list(self.starting_component.upstreams[i].name for i in range(len(self.starting_component.upstreams))))
        print(self.starting_component.inlet_fluids[0].T, self.starting_component.inlet_fluids[0].P, self.starting_component.molar_flow, self.starting_component.inlet_fluids[0].H())
        print(self.starting_component.outlet_fluids[0].T, self.starting_component.outlet_fluids[0].P, self.starting_component.molar_flow, self.starting_component.outlet_fluids[0].H())
        print("")
        
        
        
        start_index = self.components.index(self.starting_component)
        """ if you ever want to get the split count again
            #the reason we need to do this is because we can't "see ahead" with the edge_bfs method
            pre_split_count = []
            for src, node in nx.edge_bfs(G, self.components[start_index]):
                pre_split_count.append([src, node])
                
            for i in range(len(pre_split_count)):
                for j in range(len(self.components)):
                    if self.components[j] == pre_split_count[i][0]:
                        self.components[j].split_count += 1
        """
            
        
        for src, node in nx.edge_bfs(G, self.components[start_index]):
            if node.name not in visited: #or len(node.upstreams) >= 2:
                
                #this is for splitting (the component (src) before this component (node) had multiple outputs)
                for i in range(len(node.upstreams)):
                    if G.out_degree(node.upstreams[i]) >= 2 or len(node.upstreams) >= 2:
                        #print(node.upstreams[i].molar_flow, node.upstreams[i].split_count)
                        node.molar_flow += node.upstreams[i].molar_flow / G.out_degree(node.upstreams[i])
                        for j in range(len(node.upstreams[i].inlet_fluids)):
                            node.inlet_fluids.append(node.upstreams[i].outlet_fluids[j])
                
                #this is for mixing states (the current component has two inputs)
                if len(node.upstreams) >= 2:
                    mix_streams_with_upstreams(node)
                    
                #no mixing and splitting occuring
                if G.out_degree(node.upstreams[i]) < 2 and len(node.upstreams) < 2:
                    node.inlet_fluids = src.outlet_fluids
                    node.molar_flow = src.molar_flow
                    
                
                
                    
                
                
                    
                
                """
                #for mixing
                elif len(node.inlet_fluids) >= 2:
                    node.inlet_fluids = src.inlet_fluids
                    node.molar_flow = sum(src.inlet_fluids) 
                #for 
                elif len(src.upstreams) >= 2:
                    node.inlet_fluids = src.outlet_fluids
                    node.molar_flow = src.molar_flow / len(src.upstreams)
                """
                
                
                node.update_state()
                visited.append(node.name)
                
                print(node.name)
                print("current component previous connections:", list(node.upstreams[i].name for i in range(len(node.upstreams))))
                print(node.inlet_fluids[0].T, node.inlet_fluids[0].P, node.molar_flow, node.inlet_fluids[0].H())
                print(node.outlet_fluids[0].T, node.outlet_fluids[0].P, node.molar_flow, node.outlet_fluids[0].H())
                print(visited)
                print("")
                
        new_dict = []
        for src, node in nx.edge_bfs(G, self.components[start_index]):
            new_dict.append((src.name, "-->", node.name))
        print(new_dict)
    
        
            
    #def clear_results():

#wow this is really clean
#mass flow for only the input component or if you want to add more fluid (although you probably want to define the state first)
zero_compressor = Compressor(None, "compressor 0", None, 3, 0.005, 0.007, 0.65, 100000, 10000, 0, 0)

first_compressor = Compressor([zero_compressor], "compressor 1", None, 0, 0.005, 0.007, 0.65, 0, 0, 0, 10)
second_compressor = Compressor([zero_compressor], "compressor 2", None, 0, 0.005, 0.007, 0.65, 0, 100000, 0, 0)

third_compressor = Compressor([first_compressor, second_compressor, zero_compressor], "compressor 3", None, 0, 0.005, 0.007, 0.65, 50000, 0, 0, 0)
fourth_compressor = Compressor([third_compressor], "compressor 4", None, 0, 0.005, 0.007, 0.65, 50000, 0, 0, 0)
zero_compressor.upstreams = [fourth_compressor]

#print(first_compressor.enthalpy_change)

#FIX ME: using dowstreams here is problematic, while it's a lot easier to keep track of, variables cannot be accessed before they are initialized :/
#This was fixed, however I still feel like the thinking in terms of debugging while working with these is sometimes problematic, 
#however remember that this doesn't apply to how the user defines stuff, as shown above

#before_compressor_test2.downstreams = [compressor2]

constants, correlations = ChemicalConstantsPackage.from_IDs(['carbon dioxide', 'hexane'])
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)
mole_fractions = [0.5, 0.5]
state = flasher.flash(mole_fractions, T=300, P=10000)



Rankine_Cycle = System(state, flasher, zero_compressor, [zero_compressor, first_compressor, second_compressor, third_compressor, fourth_compressor])
Rankine_Cycle.run()

#print(second_compressor.inlet_fluids[0].T, second_compressor.outlet_fluids)