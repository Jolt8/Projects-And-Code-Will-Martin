import numpy as np

import math 

import scipy
from scipy import optimize

import cantera as ct

from thermo import ChemicalConstantsPackage, CEOSGas, PRMIX, CEOSLiquid, FlashVL, equilibrium, EquilibriumStream, EquilibriumState, eos_mix, Flash

from plots import plot_2d


def fenske_underwood_gilliland_method(flasher, respective_mole_fractions_in_feed, feed_molar_flow, LK_mole_fraction_in_distillate, HK_mole_fraction_in_bottoms, actual_reflux_ratio):
    #note that if actual reflux ratio is none, it will return the minimum reflux ratio so you can build on that 
    
    
    Nc = len(respective_mole_fractions_in_feed)
    
    state = flasher.flash(zs=respective_mole_fractions_in_feed, T=300, P=100000)

    
    saturation_pressures = state.Psats()
    
    component_names = state.names
    sorted_pairs = sorted(zip(saturation_pressures, component_names))
    print("EEE", sorted_pairs)
    saturation_pressures, sorted_fluids = zip(*sorted_pairs)
    saturation_pressures = list(saturation_pressures)  
    sorted_fluids = list(sorted_fluids)

    light_key_p_sat = saturation_pressures[-1] # Higher volatility
    heavy_key_p_sat = saturation_pressures[0]  # Lower volatility
    
    
    
    
    volatility_ratios = [p_sat / heavy_key_p_sat for p_sat in saturation_pressures]
    
    relative_volatility = light_key_p_sat / heavy_key_p_sat
    #usually denoted by the symmbol alpha
    
    #volatility_ratios = []
    #for i in range(0, Nc):
        #volatility_ratios.append(state.Ks(state)[i] / state.Ks(state)[-1])
        
        
    HK_mole_fraction_in_distillate = 1 - HK_mole_fraction_in_bottoms
    #print(HK_mole_fraction_in_distillate, 1 - LK_mole_fraction_in_distillate)
    LK_mole_fraction_in_bottoms = 1 - LK_mole_fraction_in_distillate
    #print(LK_mole_fraction_in_bottoms, 1 - HK_mole_fraction_in_bottoms)
    
    #Which of these two is cleaner?
    N_min = (
        math.log10((LK_mole_fraction_in_distillate / HK_mole_fraction_in_distillate) * (HK_mole_fraction_in_bottoms / LK_mole_fraction_in_bottoms))
        / (math.log10(relative_volatility))
        )

    N_min = ((math.log10((LK_mole_fraction_in_distillate / (1 - LK_mole_fraction_in_distillate)) * (HK_mole_fraction_in_bottoms / (1 - HK_mole_fraction_in_bottoms))))
                / (math.log10(relative_volatility)))
    
    #minimum_reflux_ratio 

    
    A = math.log10((1 - HK_mole_fraction_in_bottoms) / HK_mole_fraction_in_bottoms)
    #calculate the recovery ratios for all other components
    distillate_recovery_ratios = [(10 ** A) * (ratio ** N_min) / (1 + (10 ** A) * (ratio ** N_min)) for ratio in volatility_ratios]
    
    bottoms_recovery_ratios = [1 - ratio for ratio in distillate_recovery_ratios]
    
    respective_mass_flows = [comp_mole_fraction * feed_molar_flow for comp_mole_fraction in respective_mole_fractions_in_feed] 
    #(volatility_ratios[i] * distillate_recovery_ratios[i]) / (volatility_ratios[i] - underwood_constant) - 1 + state.quality
    #(volatility_ratios[i] * respective_mass_flows[i]) / (volatility_ratios[i] - optimized1.x) - 1
    
    def find_underwood_constant(underwood_constant):
        results1 = []
        for i in range(Nc):
            results1.append((volatility_ratios[i] * respective_mole_fractions_in_feed[i]) / (volatility_ratios[i] - underwood_constant))
        return sum(results1) + 1 - state.quality
    max_retries1 = 500
    retries1 = 0
    optimized1 = None
    while optimized1 is None and retries1 < max_retries1:
        optimized1 = scipy.optimize.root(find_underwood_constant, [1.5])
        retries1 += 1
    #print(optimized1)
        
    def underwood_minimum_reflux(minimum_reflux_ratio):
        results2 = []
        for i in range(Nc):
            results2.append((volatility_ratios[i] * distillate_recovery_ratios[i]) / (volatility_ratios[i] - optimized1.x))
        return sum(results2) - minimum_reflux_ratio - 1
    max_retries2 = 500
    retries2 = 0
    optimized2 = None
    while optimized2 is None and retries2 < max_retries2:
            optimized2 = scipy.optimize.root(underwood_minimum_reflux, [1.5])
            retries2 += 1
    unwrapped_value = np.ravel(optimized2.x[0])
    min_reflux_ratio = unwrapped_value[0]
    #print(optimized2)
    """
    if actual_reflux_ratio is None:
        return("compound_names", component_names, "saturation_pressures:", saturation_pressures, "relative_volatility:", relative_volatility, "volatility_ratios:", volatility_ratios, 
            "distillate_recovery_ratios:", distillate_recovery_ratios, "bottoms_recovery_ratios:", bottoms_recovery_ratios, 
            "N_min:", N_min, "underwood_constant", optimized1.x, "minimum reflux", optimized2.x,
            "light key psat", light_key_p_sat, "heavy_key psat", heavy_key_p_sat)
    """
    
    X = (actual_reflux_ratio - min_reflux_ratio) / (actual_reflux_ratio + 1)
    #print("X", X)
    
    #Y = 1 - math.exp(((1 + (54.4*X))/(11 + (117.2*X)))*((X-1)/(X)))
    #print("Y", Y)
    
    Y = (1 - (X**0.0031)) / (1 - (0.99357 * (X**0.0031)))
    print(Y)
    
    def gilland_correlation(N_actual):
        return ((N_actual - N_min) / (N_actual + 1)) - Y
    max_retries3 = 500
    retries3 = 0
    optimized3 = None
    while optimized3 is None and retries3 < max_retries3:
        optimized3 = scipy.optimize.root(gilland_correlation, [1])
        retries3 += 1
    #print(optimized3)
    unwrapped_value = np.ravel(optimized3.x[0])
    N_actual = unwrapped_value[0]
    
    if N_actual == 1:
        print("")
        print("")
        print("Your inputted actual reflux ratio:", actual_reflux_ratio, "is less than the calculated:", min_reflux_ratio, "please increase it!")
        print("")
        print("")
    
    print("compound_names", sorted_fluids, "saturation_pressures:", saturation_pressures, "volatility_ratios:", volatility_ratios,)
    print("distillate_recovery_ratios:", distillate_recovery_ratios, "bottoms_recovery_ratios:", bottoms_recovery_ratios,)      
    print("N_min:", N_min, "underwood_constant", optimized1.x, "minimum reflux", optimized2.x,)          
    print("light key psat", light_key_p_sat, "heavy_key psat", heavy_key_p_sat, )          
    print("min_reflux_ratio", min_reflux_ratio, "actual_reflux_ratio", actual_reflux_ratio, )         
    print("N_min", N_min, "N_actual", N_actual)
    print("")
    
    return("N_actual", N_actual, "min_reflux_ratio", min_reflux_ratio, "actual_reflux_ratio", actual_reflux_ratio)



constants, correlations = ChemicalConstantsPackage.from_IDs(["water", "ethanol"])
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)


print(
    fenske_underwood_gilliland_method(
        flasher=flasher,
        respective_mole_fractions_in_feed=[0.85, 0.15],
        feed_molar_flow=100,
        LK_mole_fraction_in_distillate=0.60,
        HK_mole_fraction_in_bottoms=0.60, 
        actual_reflux_ratio=1
    )
)
#Note to self, if your N_actual is 1, your inputted actual reflux ratio is probably too low

"""
returned_values = []
x_min = 0.25
x_max = 2
x_amount = 5
x = np.linspace(x_min, x_max, x_amount)

for i in range(len(x)):
    print(fenske_underwood_gilliland_method(
        flasher=flasher,
        respective_mole_fractions_in_feed=[0.4, 0.4, 0.2],
        feed_molar_flow=100,
        LK_mole_fraction_in_distillate=0.90,
        HK_mole_fraction_in_bottoms=0.90, 
        actual_reflux_ratio=x[i],
    ))
"""
"""
plot_2d(fenske_underwood_gilliland_method, 1, 5, 10, flasher, [0.6, 0.4], 100, 0.9, 0.9, None)
"""

    
    
"""
    
"""


def souders_brown_equation(state, k_value, vapor_volumetric_flow_rate, max_vapor_velocity_safety_factor):
    #note that you should input the max vapor velocity on one tray
    liq_density = state.liquid0.rho()
    gas_density = state.gas.rho()
    max_vapor_velocity = k_value * (math.sqrt((liq_density - gas_density) / gas_density))
    
    max_vapor_velocity_with_margin= max_vapor_velocity * max_vapor_velocity_safety_factor
    
    cross_sectional_area = vapor_volumetric_flow_rate / max_vapor_velocity_with_margin
    
    column_diameter = math.sqrt((4 * cross_sectional_area) / math.pi)
    
    print("Cross Sectional Area:", cross_sectional_area)
    print("Column Diameter:", column_diameter)
    
    return(column_diameter)



constants, correlations = ChemicalConstantsPackage.from_IDs(["toluene", "benzene"])
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)
state = flasher.flash(T=380.13, P=101325, zs=[0.85, 0.15])
#if the state.gas or state.liquid0 is nonetype, you need to adjust your temperature and pressure so that it's two-phase 
#this probably needs fixing because a +-1 degree difference from the actual two-phase state causes an error


def molar_flow_to_volumetric_flow(molar_flow, state):
    return molar_flow / state.rho()

#print(molar_flow_to_volumetric_flow(113, state))


souders_brown_equation(
    state=state,
    k_value = 0.107,
    vapor_volumetric_flow_rate=molar_flow_to_volumetric_flow(2.5, state),
    max_vapor_velocity_safety_factor=0.8
)
