import CoolProp.CoolProp as cp
from CoolProp.CoolProp import PropsSI
import scipy.optimize

import ht as ht 
from ht import insulation, conduction, conv_external, conv_internal, boiling_plate

#from CoolProp.CoolProp import IProps, get_Fluid_in dex
#from CoolProp.CoolProp import PhaseSI
from CoolProp.HumidAirProp import HAPropsSI

import colebrook

import numpy as np

import math 

import scipy
from scipy import optimize

from pyfluids import Fluid, FluidsList, Input 
from plots  import plot_3d, plot_2d

def gibbs_free_energy(heat_of_reaction, temperature, entropy_change):
    return heat_of_reaction - (temperature * entropy_change)
print(gibbs_free_energy(-92000, 650, -200))


def vant_hoft(heat_of_reaction, temperature, entropy_change):
    ideal_gas_constant = 8314
    return np.exp(-1 * (heat_of_reaction / (temperature * ideal_gas_constant))) * np.exp(entropy_change / ideal_gas_constant)

print(vant_hoft(-92000, 723, -200))
plot_2d(vant_hoft, 1, 100000, 1000, 0.0000001, 1, -92000, None, -200)
    
def Erying(heat_of_reaction, delta_enthalpy, temeperature):
    #temeperature is in kelvin
    ideal_gas_constant = 8.314
    return np.exp((-(heat_of_reaction / (ideal_gas_constant * temeperature)) + (delta_enthalpy / ideal_gas_constant)))

print("Erying", Erying(-92000, -200, 723))
#plot_3d(Erying, -92000, 92000, 800, -200, 200, 800, None, None, 723)


def arrhenius(pre_exponential_factor, activation_energy, temperature):
    ideal_gas_constant = 8.314
    return pre_exponential_factor * np.exp(activation_energy / temperature * ideal_gas_constant)

plot_3d(arrhenius, -1000, 100000, 200, 0, 1000, 200, None, 460, None)
            
    #-200, 200, 200, 0, 600, None, 460, None)

    
    