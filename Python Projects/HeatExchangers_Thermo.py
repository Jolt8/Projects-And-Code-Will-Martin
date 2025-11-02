import ht 

import numpy as np

import math 

import scipy
from scipy import optimize

from thermo import ChemicalConstantsPackage, CEOSGas, PRMIX, CEOSLiquid, FlashVL, equilibrium, EquilibriumStream, EquilibriumState, eos_mix, Flash



    
    
constants, correlations = ChemicalConstantsPackage.from_IDs(["water", "hexane"])
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)
hot_state = flasher.flash(T=370, P=100000, zs=[0, 1])
cold_state = flasher.flash(T=300, P=100000, zs=[1, 0])


print(effectiveness_NTU_thermo(hot_stream_state=hot_state, hot_stream_mass_flow=10, cold_stream_state=cold_state, cold_stream_mass_flow=20,
                               subtype="counterflow", hot_inlet_temp=370, hot_outlet_temp=350, cold_inlet_temp=300, cold_outlet_temp=None, UA=None, n_shell_tube=None))