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
def H_dot_FTPQ(*args):
    def pre_H_dot_FTPQ(fluid, H_dot, T1, T2, P1, P2, Q1, Q2):
        """Function that Calculates a change in Enthalpy 
        

        Args:
            fluid (fluid): Any fluid in the coolprop fluid list
            H_dot (Kj/kg): Enthalpy Change
            T1 (K): starting temperature
            T2 (K): ending temperature
            P1 (p): startng pressure
            P2 (p): ending pressure
            Q1 (0 to 1): starting vapour quality
            Q2 (0 to 1): ending vapour quality

        Returns:
            Any input with "None": 
        """
        if T1 is None and T2 is None:
            return (cp.PropsSI("H", "P", P2, "Q", Q2, fluid) - cp.PropsSI("H", "P", P1, "Q", Q1, fluid) - H_dot)
        elif P1 is None and P2 is None:
            return cp.PropsSI("H", "T", T2, "Q", Q2, fluid) - cp.PropsSI("H", "T", T1, "Q", Q1, fluid) - H_dot
        else: 
            return cp.PropsSI("H", "T", T2 + 400, "P", P2, fluid) - cp.PropsSI("H", "T", T1, "P", P1, fluid) - H_dot
    #solve it
    i = args.index(None)
    return scipy.optimize.fsolve(lambda x: pre_H_dot_FTPQ(*args[:i], x, *args[i+1:]), 1)

##print("Enthalpy Change", H_dot_FTPQ("Water", 100000, 300, None, 100000, 100000, 0, 1))

# Latent heat function
def latent_heat(fluid):
    return (cp.PropsSI("H", "T", 300, "Q", 1, fluid) - cp.PropsSI("H", "T", 300, "Q", 0, fluid))

##print("Latent heat:", latent_heat("Water"))


#Mass Flow Calculations

def mass_flow_FT1T2W(*args):
    def pre_mass_flow_FWT(m_dot, fluid, wattage, T_in, T_out):
        """
        Gets mass flow rate for required wattage based on input and output temperature
        Args:
            m_dot (kg/s): _description_
            fluid (): Any fluid from the coolprop fluid list
            T_in (K): temperature of fluid going in
            T_out (K): temperature of fluid going out
            W (W): required wattage
            
        Intermediates:
            delta_T (K): Change in temperature
            
        Returns:
           Any input with "None":
        """
        
        delta_T = T_out - T_in
        average_C = (cp.PropsSI("C", "T", T_out, "Q", 0, fluid) + cp.PropsSI("C", "T", T_in, "Q", 0, fluid)) / 2
        return (m_dot * (average_C * delta_T) - wattage)
    i = args.index(None)
    return scipy.optimize.fsolve(lambda x: pre_mass_flow_FWT(*args[:i], x, *args[i+1:]), 1)
    
##print(mass_flow_FT1T2W(0.26, "Water", None, 300, 320))


assumed_final_temperature = 85 + 273.15
def T_final_F1TmPQ_F2TmPQ(*args):
    def pre_final_T_FTmPQ_FTmPQ(final_T, fluid1, T1, m1, P1, Q1, fluid2, T2, m2, P2, Q2):
        """Finds the final tempearture even if the fluid boils 

        Args:
            final_T (K): final temperature of both fluids
        
            fluid1 (str): fluid 1 ex. "Water"
            T1 (K): temperature
            m1 (kg or kg/s): mass or mass flow rate
            P1 (p): pressure
            Q1 (0-1): vapor quality (0 = sat liquid, 1 = sat vapour)
            
            fluid2 (str): fluid 2 ex. "Water"
            T2 (K): temperature 
            m2 (kg or kg/s): mass or mass flow rate 
            P2 (p): pressure
            Q2 (0-1): vapor quality (0 = sat liquid, 1 = sat vapour)
            
        intermediates:
            cp (kj/kg/K): specific heat
            T_final_no_cp (): final temperature without latent heat
            bp (p): boiling point of the fluid at a specific pressure
            h_vap (p): latent heat of vaporization 
        
        Returns:
            Any value except for fluid provided with "None"
        """
        cp1 = cp.PropsSI("C", "T", T1, "Q", Q1, fluid1)  #specific heat 1
        cp2 = cp.PropsSI("C", "T", T2, "Q", Q2, fluid2)  #specific heat 2
        
        T_final_no_cp = (m1 * cp1 * T1 + m2 * cp2 * T2) / (m1 * cp1 + m2 * cp2)  #final temperature without latent heat 
        
        bp1 = cp.PropsSI("T", "P", P1, "Q", Q1, fluid1)  #boiling point at pressure and quality of fluid 1
        bp2 = cp.PropsSI("T", "P", P2, "Q", Q2, fluid2)  #boiling point at pressure and quality of fluid 2 

        h_vap1 = (cp.PropsSI("H", "T", T1, "Q", 1, fluid1) - cp.PropsSI("H", "T", T1, "Q", Q1, fluid1))  #latent heat 1
        h_vap2 = (cp.PropsSI("H", "T", T2, "Q", 1, fluid2) - cp.PropsSI("H", "T", T2, "Q", Q1, fluid2))  #latent heat 2

        if bp1 >= T_final_no_cp and bp2 >= T_final_no_cp:  #if both fluids don't boil
            return ((m1 * cp1 * T1 + m2 * cp2 * T2) / (m1 * cp1 + m2 * cp2) - final_T)
        elif bp1 < T_final_no_cp and bp2 < T_final_no_cp:  #if both fluids boil
            return ((m1 * h_vap1 * T1 + m2 * h_vap2 * T2) / (m1 * h_vap1 + m2 * h_vap2) - final_T)
        elif bp1 <= T_final_no_cp:  #if fluid1 boils
            return ((m1 * h_vap1 * T1 + m2 * cp2 * T2) / (m1 * h_vap1 + m2 * cp2) - final_T)
        else:  #if fluid2 boils
            return ((m1 * cp1 * T1 + m2 * h_vap2 * T2) / (m1 * cp1 + m2 * h_vap2) - final_T)
    i = args.index(None)
    return scipy.optimize.fsolve(lambda x: pre_final_T_FTmPQ_FTmPQ(*args[:i], x, *args[i+1:]), 300)

#print("solvetest", T_final_F1TmPQ_F2TmPQ (300, "Water", 368, 0.1, 100000, 0, "Isobutane", None, 0.1, 1200000, 0))


def find_missing_properties (fluid, T, P, Q, rho):
#rho (density) is refered to as D in coolprop 
    properties = [T, P, Q, rho]
    property_names = ["T", "P", "Q", "D"]
    
    missing_properties_positions = []
    missing_property_names = []
    given_properties_positions = []
    given_property_names = []
    
    for i, prop in enumerate(properties):
        if prop is None:
            missing_properties_positions.append(i)
            missing_property_names.append(property_names[i])
        else:
            given_properties_positions.append(i)
            given_property_names.append(property_names[i])
            
    if len(missing_properties_positions) == 0:
        raise ValueError("No properties are missing; nothing to calculate.")
    if len(missing_properties_positions) > 2:
        raise ValueError("Too many missing properties; only one or two can be calculated.")
    
    
    for i, missing_index in enumerate(missing_properties_positions):
        #These if and elif statements deal with the issue of getting Q (Phase/Quality) from just T and P or D and P
        if given_properties_positions == [0, 1] and 3 in missing_properties_positions:
            properties[2] = cp.PhaseSI("T", T, "P", P, fluid)
            properties[3] = cp.PropsSI("D", "T", T, "P", P, fluid)
            continue
        elif given_properties_positions == [1, 3] and 3 in missing_properties_positions:
            cp.PhaseSI("P", P, "Q", Q, fluid)
            properties[2] = cp.PhaseSI("P", P, "D", rho, fluid)
            properties[3] = cp.PropsSI("P", "T", T, "D", rho, fluid)
            continue
        
        # Extract the CoolProp-compatible names and values
        input_name_1 = given_property_names[0]
        input_value_1 = properties[given_properties_positions[0]]
        
        input_name_2 = given_property_names[1]
        input_value_2 = properties[given_properties_positions[1]]

        # Calculate the missing property
        calculated_value = cp.PropsSI(
            property_names[missing_index],
            input_name_1, input_value_1,
            input_name_2, input_value_2,
            fluid
        )
        
    #Coolprop Can't get Q values with only T and P or D and P 
        
        # Update the properties list with the calculated value
        properties[missing_index] = calculated_value
    if properties[2] == -1.0:
            raise "returned a Q value of -1, consider rasing your temperature, lowering your pressure, or raising your density"
    return properties 
#print (find_missing_properties("Water", None, 1000, None, 0.1))
#print (find_missing_properties("Isobutane", 300, 100000, None, None))
#print (find_missing_properties("Water", None, 100000, None, 920))

def find_Q_final_with_Tsat(fluid, saturation_temperature, actual_pressure):
    """
    Keep in mind that to find quality you must specify a position on the saturation curve and then use an increase in entropy or pressure to find the quality 
    """
    saturation_pressure = cp.PropsSI("P", "T", saturation_temperature, "Q", 0, fluid)
    Hliq = cp.PropsSI("H", "T", saturation_temperature, "Q", 0, fluid)
    
    #print (saturation_pressure)
    #find the difference between the the actual pressure and the saturation pressure
    Hinput = cp.PropsSI("H", "Q", 0, "P", actual_pressure, fluid) - cp.PropsSI("H", "Q", 0, "P", saturation_pressure, fluid) 
    H = Hliq + Hinput
    #print(cp.PropsSI("H", "Q", 0, "P", actual_pressure, fluid), Hliq, Hinput, H)
    Q = cp.PropsSI("Q", "P", saturation_pressure, "H", H, fluid)
    
    if Q == -1 and actual_pressure < saturation_pressure:
        return 1
    elif Q == -1 and actual_pressure >= saturation_pressure:
        return 0
    else:
        return Q

#print("function final Q Tsat", find_Q_final_with_Tsat("Isobutane", 20 + 273.13, 1500000))

def find_Q_final_with_Psat(fluid, saturation_pressure, actual_pressure):
    """
    Keep in mind that to find quality you must specify a position on the saturation curve and then use an increase in entropy or pressure to find the quality 
    """
    Hliq = cp.PropsSI("H", "P", saturation_pressure, "Q", 0, fluid)
    
    #find the difference between the the actual pressure and the saturation pressure
    Hinput = cp.PropsSI("H", "Q", 0, "P", actual_pressure, fluid) - cp.PropsSI("H", "Q", 0, "P", saturation_pressure, fluid) 
    H = Hliq + Hinput
    #print(cp.PropsSI("H", "Q", 0, "P", actual_pressure, fluid), Hliq, Hinput, H)
    Q = cp.PropsSI("Q", "P", saturation_pressure, "H", H, fluid)
    
    
    if Q == -1 and actual_pressure < saturation_pressure:
        return 1
    elif Q == -1 and actual_pressure >= saturation_pressure:
        return 0
    else:
        return Q
    return Q

#print("function final Q Psat", find_Q_final_with_Psat("Isobutane", 300000, 1500000))


def find_Q_final_with_dH(fluid, input_enthalpy, output_enthalpy, saturation_pressure):
    """
    Keep in mind that to find quality you must specify a position on the saturation curve and then use an increase in entropy or pressure to find the quality 
    """
    
    H = output_enthalpy
    #print(cp.PropsSI("H", "P", saturation_pressure, "Q", 0, fluid))
    Q = cp.PropsSI("Q", "P", saturation_pressure, "H", H, fluid)
    return Q

#print("function final Q dH", find_Q_final_with_dH("Isobutane", 300000, 304000, 300000))

#print(cp.PropsSI("P", "Q", 1, "T", 293, "isobutane"))

