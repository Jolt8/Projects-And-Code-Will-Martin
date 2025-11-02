from CoolProp.CoolProp import PropsSI, get_phase_index, PhaseSI
import CoolProp.CoolProp as cp
from CoolProp.HumidAirProp import HAPropsSI

import numpy as np

import math 

import scipy
from scipy import optimize

from operator import itemgetter

from pyfluids import Fluid, FluidsList, Input 
def gaseous_fluid_list_T1_T2_P1_P2(temperature, pressure):
    """Finds all fluids that are gaseous at a temperature and pressure

    Args:
        temperature (K): lowest acceptible temperature
        pressure (P): highest acceptible pressure
        Remember that higher pressure = less gas, more liquid

    Returns:
        gaseous_fluids (list): list of gaseous fluids
    """
    gaseous_fluids = []
    Fluids = list(cp.FluidsList())
    for i in range(len(Fluids)):
        try:
            if cp.PhaseSI("T", temperature, "P", pressure, list(cp.FluidsList())[i]) == "liquid":
                gaseous_fluids.append(Fluids[i])
        except:
            pass
        i = i + 1
    return gaseous_fluids
#print("gaseous fluids", gaseous_fluid_list_T1_T2_P1_P2(300, 100000))


def sweet_spot_fluids_from_range(lowtemperature, hightemperature, lowpressure, highpressure):
    """Finds all fluids that are liquid at temperature2 and pressure2 and gaseous at temperature2 and pressure2

    Args:
        temperature1 (K): lowest accepted temperature
        temperature2 (K): highest accepted temperature
        pressure1 (p): lowest accepted pressure
        pressure2 (p): highest accepted pressure
        
    Intermediates:
        liquid_test (str): outputs "gas" if something is a gas
        gas_test (str): outputs "liquid" if something is a liquid

    Returns:
        sweet_fluids (list): list of fluids that are gaseous and liquid within this temperature range
    """
    sweet_fluids = []
    Fluids = list(cp.FluidsList())
    for i in range(len(Fluids)):
        try:
            liquid_test = cp.PhaseSI("T", lowtemperature, "P", highpressure, list(cp.FluidsList())[i])
            gas_test = cp.PhaseSI("T", hightemperature, "P", lowpressure, list(cp.FluidsList())[i])
            if liquid_test == "liquid" and gas_test == "gas":
                fluid_with_props = [Fluids[i]]
                sweet_fluids.append(fluid_with_props)
        except:
            pass
    return sweet_fluids
#print("sweet_spot_range", sweet_spot_fluids_from_range(300, 350 , 100000, 400000))



""" 
try: boiling_point_1 = cp.PropsSI("T", "P", lowpressure, "Q", 0, list(cp.FluidsList())[i]) 
        except: pass
        try: boiling_point_2 = cp.PropsSI("T", "P", highpressure, "Q", 0, list(cp.FluidsList())[i])
        except: pass
        try: vapour_pressure_1 = cp.PropsSI("P", "T", lowtemperature, "Q", 0, list(cp.FluidsList())[i])
        except: pass
        try: vapour_pressure_2 = cp.PropsSI("P", "T", hightemperature, "Q", 0, list(cp.FluidsList())[i])
        except: pass

"""

def sweet_spot_fluids_optimized(ideal_temperature, above_T, below_T, ideal_pressure, above_P, below_P):
    sweet_fluids = []
    Fluids = list(cp.FluidsList())
    for i in range (len(Fluids)):
        # Try to find the boiling point
        try: boiling_point = cp.PropsSI("T", "P", ideal_pressure, "Q", 0, Fluids[i]) 
        except: pass
        
        # Try to find the vapour pressure
        try: vapour_pressure = cp.PropsSI("P", "T", ideal_temperature, "Q", 0, Fluids[i])
        except: pass
        
        # See if the fluid is a liquid at the lowest temperature and highest pressure
        # See if the fluid is a gas at the highest temperature and lowest pressure
        try:
            liquid_test = cp.PhaseSI("T", ideal_temperature - below_T, "P", ideal_pressure + above_P, Fluids[i])
            gas_test = cp.PhaseSI("T", ideal_temperature + above_T, "P", ideal_pressure - below_P, Fluids[i])
            
            # If it can exist as a liquid and a gas at the given temperatures and pressures, put it in the list
            if liquid_test == "liquid" and gas_test == "gas":
                fluid_with_props = [Fluids[i], boiling_point, vapour_pressure]
                sweet_fluids.append(fluid_with_props)
        except:
            pass
        i = i + 1
    sort_T_asc = sorted(sweet_fluids, key=lambda x: abs(x[1] - ideal_temperature))
    sort_P_asc = sorted(sweet_fluids, key=lambda x: abs(x[2] - ideal_pressure))
    new_list = []
    for i in range (len(sort_T_asc)):
        new_list.append(sort_T_asc[i][0])
    return new_list
#print ("sweet_spot_optimized", sweet_spot_fluids_optimized(300, 10, 10, 300000, 150000, 150000))
#print ("sweet_spot_optimized", sweet_spot_fluids_optimized(380, 10, 10, 1500000, 150000, 150000))

def property_analysis_universal(property, type="closest", ideal_property_value=100, min=0, max=100000, temperature=300, pressure=100000):
    sweet_fluids = []
    Fluids = list(cp.FluidsList())
   
    for i in range (len(Fluids)):
        try: 
            if cp.PropsSI(property, "T", temperature, "P", pressure, Fluids[i]) >= min and cp.PropsSI(property, "T", temperature, "P", pressure, Fluids[i]) <= max:
                sweet_fluids.append([Fluids[i], cp.PropsSI(property, "T", temperature, "P", pressure, Fluids[i])])
        except:
            sweet_fluids.append([Fluids[i], 0])
    if type == "closest":
        return sorted(sweet_fluids, key=lambda x:abs(ideal_property_value-x[1]))
    elif type == "min":
        return sorted(sweet_fluids, key=lambda x: x[1], reverse=False)
    elif type == "max":
        return sorted(sweet_fluids, key=lambda x: x[1], reverse=True)
    else:
        raise "Invalid type, use closest, min, or max"
    



def trivial_property_analysis_closest(property, ideal_property_value, allowed_above, allowed_below):
    """Returns a list of fluids that are listed based on how close they are to the ideal value
    

    Args:
        property (): coolprop compatible property
        order (how the list should be ordered): 0 = closest to idea, 1 = least, 2, = greatest 
        ideal_property_value (_type_): _description_
        allowed_above (_type_): _description_
        allowed_below (_type_): _description_

    Returns:
        _type_: _description_
    """
    sweet_fluids = []
    Fluids = list(cp.FluidsList())
    min_value = ideal_property_value - allowed_below
    max_value = ideal_property_value + allowed_above
   
    for i in range (len(Fluids)):
        try: 
            if cp.PropsSI(property, Fluids[i]) >= min_value and cp.PropsSI(property, Fluids[i]) <= max_value:
                sweet_fluids.append([Fluids[i], cp.PropsSI(property, Fluids[i])])
        except:
            sweet_fluids.append([Fluids[i], 0])
    return sorted(sweet_fluids, key=lambda x:abs(ideal_property_value-x[1]))
#print("closest trivial property analysis", trivial_property_analysis_closest("GWP100", 100, 10000, 10000))

def trivial_property_analysis_least(property, min_value, max_value):
    sweet_fluids = []
    Fluids = list(cp.FluidsList())
    for i in range (len(Fluids)):
        try: 
            if cp.PropsSI(property, Fluids[i]) >= min_value and cp.PropsSI(property, Fluids[i]) <= max_value:
                sweet_fluids.append([Fluids[i], cp.PropsSI(property, Fluids[i])])
        except:
            sweet_fluids.append([Fluids[i], 0])
    return sorted(sweet_fluids, key=lambda x: x[1], reverse=False)
#print("least trivial property analysis", trivial_property_analysis_least("GWP100", 0, 100))

def trivial_property_analysis_greatest(property, min_value, max_value):
    sweet_fluids = []
    Fluids = list(cp.FluidsList())
    for i in range (len(Fluids)):
        try: 
            if cp.PropsSI(property, Fluids[i]) >= min_value and cp.PropsSI(property, Fluids[i]) <= max_value:
                sweet_fluids.append([Fluids[i], cp.PropsSI(property, Fluids[i])])
        except:
            sweet_fluids.append([Fluids[i], 0])
    return sorted(sweet_fluids, key=lambda x: x[1], reverse=True)
#print("greatest trivial property analysis", trivial_property_analysis_greatest("GWP100", 0, 100))

    



def choose_best_fluid(evaporator_list, condenser_list):
    # Create dictionaries to store the ranks of each fluid in both lists
    evaporator_rank = {fluid: rank + 1 for rank, fluid in enumerate(evaporator_list)}
    condenser_rank = {fluid: rank + 1 for rank, fluid in enumerate(condenser_list)}

    # Combine the ranks: calculate a combined score for each fluid
    possible_fluids = (np.unique(evaporator_list + condenser_list))
    best_fluids = []
    for fluid in possible_fluids:
        evap_rank = evaporator_rank.get(fluid, float("inf"))  # Use a high rank if fluid is missing
        cond_rank = condenser_rank.get(fluid, float("inf"))  # Same for condenser list
        best_fluids.append([fluid, evap_rank, cond_rank])
        

    best_fluids = sorted(best_fluids, key=lambda x: (x[1] + x[2]))
    return best_fluids
    
    
    sort_key = {}
    best_fluids = {}
    for fluid in evaporator_list:
        
        sort_key[fluid] = evap_rank + cond_rank
        best_fluids[fluid] = evap_rank, cond_rank
    
    # Sort fluids by combined rank
    sorted(sort_key.items(), key=lambda x: x[1])
    return best_fluids

#print("best fluid", choose_best_fluid(sweet_spot_fluids_optimized(90 + 273, 20, 20, 1500000, 300000, 300000), sweet_spot_fluids_optimized(20 + 273, 20, 20, 300000, 150000, 150000)))

print(cp.PropsSI("H", "T", 300, "P", 100000, "water"))
print(cp.PropsSI("H", "T", 300, "P", 300000, "water"))

 
def fluid_optimizer(evap_weight, 
                    evap_temp, evap_above_T, evap_below_T, 
                    evap_pressure, evap_above_P, evap_below_P,
                    cond_weight,
                    cond_temp, cond_above_T, cond_below_T, 
                    cond_pressure, cond_above_P, cond_below_P,
                    properties=[["C", 1, "max", 4000, 0, float("inf")], ["GWP100", 0.5, "min", 400, 0, 3000]]
                    ):
    """The properties are formatted like this [coolprop property name, property weight, type ]
    """
    # Get potential fluids for evaporator and condenser
    evap_list = sweet_spot_fluids_optimized(evap_temp, evap_above_T, evap_below_T, evap_pressure, evap_above_P, evap_below_P)
    cond_list = sweet_spot_fluids_optimized(cond_temp, cond_above_T, cond_below_T, cond_pressure, cond_above_P, cond_below_P)
    
    # Create ranking dictionaries
    evaporator_rank = {fluid: rank + 1 for rank, fluid in enumerate(evap_list)}
    condenser_rank = {fluid: rank + 1 for rank, fluid in enumerate(cond_list)}
    
    temperature = (evap_temp + cond_temp) / 2
    pressure = (evap_pressure + cond_pressure) / 2
    
    # Analyze properties
    property_lists = [property_analysis_universal(prop[0], prop[2], prop[3], prop[4], prop[5], temperature, pressure) for prop in properties]
    
    # Extract fluid names for properties
    property_ranks = []
    for prop_list in property_lists:
        ranked_fluids = {fluid[0]: rank + 1 for rank, fluid in enumerate(prop_list)}
        property_ranks.append(ranked_fluids)
    
    # Calculate total weight factors
    other_property_total_weight = sum(prop[1] for prop in properties)
    total_weight = evap_weight + cond_weight + other_property_total_weight
    
    evap_factor = evap_weight / total_weight
    cond_factor = cond_weight / total_weight
    property_factors = [prop[1] / total_weight for prop in properties]
    
    # Get unique fluids
    if properties is not None:
        possible_fluids = np.unique(evap_list + cond_list + [fluid[0] for prop_list in property_lists for fluid in prop_list])
    if properties is not None:
        possible_fluids = np.unique(evap_list + cond_list + [fluid[0] for prop_list in property_lists for fluid in prop_list])
    
    best_fluids = []
    for fluid in possible_fluids:
        evap_rank = evaporator_rank.get(fluid, float("inf"))
        cond_rank = condenser_rank.get(fluid, float("inf"))
        
        fluid_property_ranks = [prop_rank.get(fluid, float("inf")) for prop_rank in property_ranks]
        
        best_fluids.append([fluid, evap_rank, cond_rank] + fluid_property_ranks)
    
    # Sort based on weighted ranks
    best_fluids.sort(key=lambda x: (evap_factor * x[1] + cond_factor * x[2] + sum(pf * x[i+3] for i, pf in enumerate(property_factors))))
    
    return best_fluids

# Example usage
print("final_fluid_optimizer",
    fluid_optimizer(
        evap_weight=10,
        evap_temp=333, 
        evap_above_T=20, 
        evap_below_T=20, 
        evap_pressure=2000000, 
        evap_above_P=100000, 
        evap_below_P=100000,
        
        cond_weight=10,
        cond_temp=293, 
        cond_above_T=20, 
        cond_below_T=20, 
        cond_pressure=750000, 
        cond_above_P=100000, 
        cond_below_P=100000,
        properties=[["C", 0, "max", 4000, 0, 10000]]#, ["GWP100", 0, "min", 0, 0, 10000], ["FH", 0, "min", 0, 0, 10000]]
    )
)