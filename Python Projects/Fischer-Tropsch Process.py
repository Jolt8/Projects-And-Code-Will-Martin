import CoolProp.CoolProp as cp
from CoolProp.CoolProp import PropsSI
import scipy.odr
import scipy.optimize

#from PressureDropCalculations import reynolds_number
from PressureDropCalculations import reynolds_number

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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

from plots import plot_3d

def final_solver1(chain_growth, CO_percentage, H2_percentage, temperature):
        return (0.2332 * (CO_percentage / (H2_percentage))) + 0.633 * (1 - 0.0039 *((temperature + 273) - 533)) - chain_growth
    
    
def final_solver2(CO_percentage, H2_percentage, temperature):
    return (0.2332 * (CO_percentage / (H2_percentage))) + 0.633 * (1 - 0.0039 *((temperature + 273) - 533))

def final_solver3(CO_percentage, temperature):
    return (0.2332 * (CO_percentage)) + 0.633 * (1 - 0.0039 *((temperature + 273) - 533))

plot_3d(final_solver3, 0, 1.5, 50, 200, 500, 50, None, None)


def test(*args):
    i = args.index(None)
    return scipy.optimize.fsolve(lambda x: final_solver1(*args[:i], x, *args[i+1:]), 1.5)
print("test", test(1, 1, 1, None))

def fischer_hydrocarbon_ratios(chain_growth, target_carbon_number):
    carbon_list = []
    chain_growth_float = chain_growth
    for i in range(40):
        mass_fraction = (chain_growth_float ** (i - 1)) * ((1 - chain_growth_float) ** 2) * i
        carbon_list.append(["C", i, mass_fraction])
    return carbon_list
    #sorted(carbon_list, key=lambda x:abs(target_carbon_number-x[2]))
print("test", fischer_hydrocarbon_ratios(0.8, 10))

def fischer_hydrocarbon_range (min_hydrocarbon, max_hydrocarbon):
    ratio_list = fischer_hydrocarbon_ratios(0.8, 10)
    concentration_total = 0
    for i in range(len(ratio_list)):
        if ratio_list[i][1] >= min_hydrocarbon and ratio_list[i][1] <= max_hydrocarbon:
            concentration_total = concentration_total + ratio_list[i][2]
    return concentration_total
print("hydrocarbon_range", fischer_hydrocarbon_range(6, 12))

def fischer_hydrocaron_ratios_optimized(chain_growth, target_carbon_number, mass_fraction):
    return ((chain_growth ** (target_carbon_number - 1)) * ((1 - chain_growth) ** 2) * target_carbon_number) - mass_fraction

def test2(*args):
    i = args.index(None)
    return scipy.optimize.fsolve(lambda x: fischer_hydrocaron_ratios_optimized(*args[:i], x, *args[i+1:]), 0.5)
print("Final test", test2(None, 9, 0.06))

