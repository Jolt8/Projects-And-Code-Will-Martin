import numpy as np

import math 

import colebrook

import ht as ht 

import scipy
from scipy import optimize

import cantera as ct

from thermo import ChemicalConstantsPackage, CEOSGas, PRMIX, CEOSLiquid, FlashVL, equilibrium, EquilibriumStream, EquilibriumState, eos_mix, Flash


constants, correlations = ChemicalConstantsPackage.from_IDs(['Water', 'Mercury'])
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)
mole_fractions = [0.5, 0.5]


def reynolds_number(state, fluid_velocity, characteristic_length):
    fluid_density = state.bulk.rho()
    dynamic_viscosity = state.bulk.mu()
    return (fluid_density * fluid_velocity * characteristic_length) / (dynamic_viscosity)

def grashof_number(state, characteristic_length, heat_exchange_surface_temperature):
    #I think this might be the main cause of low values
    gravitational_acceleration = 9.81
    coefficient_of_expansion = state.isobaric_expansion()
    
    dynamic_viscosity = state.mu()
    density = state.rho()
    kinematic_viscosity = dynamic_viscosity / density
    
    #the abs is to make sure it works for cases where wall temperature are lesser than ambient fluid temperature
    return (gravitational_acceleration * coefficient_of_expansion * (abs(heat_exchange_surface_temperature - state.T)) * (characteristic_length ** 3)) / (kinematic_viscosity ** 2)


def ht_internal_conv(state, state_mass_flow, pipe_characteristic_length, pipe_length, absolute_pipe_roughness):
    """Finds the nusselt number for internal convection (convection from liquid to the pipe)"""
    
    fluid_velocity = pipe_characteristic_length / (state.bulk.rho() * state_mass_flow)
    Re = reynolds_number(state, fluid_velocity, pipe_characteristic_length)
    prandtl = state.Prandtl()
    relative_pipe_roughness = (absolute_pipe_roughness / pipe_characteristic_length)
    friction_factor = colebrook.bntFriction(Re, relative_pipe_roughness, 4)
    return ht.conv_internal.Nu_conv_internal(Re, prandtl, relative_pipe_roughness, pipe_characteristic_length, pipe_length, friction_factor)
#print ("pipe_heat_loss_internal", ht_internal_conv("Water", 300, 100000, 3, 0.0125, 3, 0.0013))


def ht_external_pipe(state, external_state, external_fluid_temp, external_fluid_pressure, fluid_velocity, pipe_diameter):
    """Finds the nusselt number for external convection (convection from the pipe to a medium outside of the pipe)"""
    Re = reynolds_number(state, external_fluid_temp, external_fluid_pressure, fluid_velocity, pipe_diameter)
    
    fluid_prandtl = state.Prandtl()
    
    external_state = flasher.flash(zs=state.zs, T=external_fluid_temp, P=external_fluid_pressure)
    external_prandtl = external_state.Prandtl()
    
    fluid_dynamic_viscosity = state.mu()
    #TODO: investigate these two similar values
    external_dynamic_viscosity = state.mu()
    return ht.conv_external.Nu_external_cylinder(Re, fluid_prandtl, external_prandtl, fluid_dynamic_viscosity, external_dynamic_viscosity)
#print ("pipe_heat_loss_external", ht_external_pipe("Water", 300, 100000, 3, 0.01))


def ht_external_plate(state, fluid_velocity, pipe_characteristic_length):
    """Finds the nusselt number for internal convection (convection from plate heat exchanger to a medium outside of pipe)"""
    
    Re = reynolds_number(state, fluid_velocity, pipe_characteristic_length)
    fluid_prandtl = state.Prandtl()
    return ht.conv_external.Nu_external_horizontal_plate(Re, fluid_prandtl, pipe_characteristic_length)
#print("plate_heat_loss_external", ht_external_plate("Water", 300, 100000, 3, 0.01))

def ht_submerged_coil(ambient_state, heat_exchange_surface_temperature, outer_characteristic_length, horizontal=False):
    """Finds the nusselt number for pipes submerged in a fluid"""
    Gr = grashof_number(ambient_state, outer_characteristic_length, heat_exchange_surface_temperature)
    fluid_prandtl = ambient_state.Prandtl()
    return ht.conv_free_immersed.Nu_coil_Xin_Ebadian(fluid_prandtl, Gr, horizontal)
#print("Submerged coil heat loss external", ht_submerged_coil("Water", 350, 100000, 400, 0.125, False))


def ht_boiling_flow_chen_bennet(state, mass_flow_rate, 
                                fluid_temp_at_wall_excess, fluid_pressure_at_wall_excess,
                                pipe_characteristic_length, 
                                wall_temp
                                ):
    try:
        fluid_phase_density = state.liquid0.rho()
        gas_density = state.gas.rho()
    
    except Exception as e:
        print("Fluid in not two-phase! State Phase:", state.phase, e)
    
    fluid_phase_viscosity = state.liquid0.mu()
    gas_phase_viscosity = state.gas.mu()
    
    fluid_phase_thermal_conductivity = state.liquid0.k()
    fluid_phase_specific_heat_capacity = state.liquid0.Cp()
    fluid_phase_heat_of_vaporization = sum(wall_temp_state.Hvaps()) / (len(wall_temp_state.Hvaps()))
    fluid_surface_tension = state.sigma()
    
    quality_at_exchanger_specific_interval = state.quality
    
    wall_temp_state = flasher.flash(T=fluid_temp_at_wall_excess, P=fluid_pressure_at_wall_excess, zs=state.zs)
    
    state_avg_psat = sum(state.Psats()) / (len(state.Psats()))
    wall_state_avg_psat = sum(wall_temp_state.Psats()) / (len(wall_temp_state.Psats()))

    
    delta_saturation_pressure =  wall_state_avg_psat - state_avg_psat
    
    
    """
    volumetric_flow_rate = fluid_velocity * ((pipe_characteristic_length ** 2) * math.pi)
    mass_flow_rate = volumetric_flow_rate * fluid_phase_temperature
    """
    
    return ht.boiling_flow.Chen_Bennett(mass_flow_rate, quality_at_exchanger_specific_interval, pipe_characteristic_length, 
                                        fluid_phase_density, gas_density,
                                        fluid_phase_viscosity, gas_phase_viscosity,
                                        fluid_phase_thermal_conductivity, fluid_phase_specific_heat_capacity, fluid_phase_heat_of_vaporization,
                                        fluid_surface_tension, delta_saturation_pressure, 
                                        wall_temp
                                        )


def heat_transfer_coeff(state, Nu, characteristic_length):
    fluid_thermal_cond = state.bulk.k()
    return (Nu * fluid_thermal_cond) / characteristic_length

def overall_heat_transfer_coeff(internal_nusselt_number, internal_state, internal_characteristic_length, internal_area, 
                                external_nusselt_number, external_state, external_characteristic_length, external_area, 
                                layers):
    """
    Keep in mind that layers if formatted [(thickness, thermal conductivity, area), (thickness, thermal conductivity, area)]
    """
    conductive_resistance = sum(thickness / (conductivity * area) for thickness, conductivity, area in layers)
    
    internal_heat_transfer_coeff = heat_transfer_coeff(internal_state, internal_nusselt_number, internal_characteristic_length)
    external_heat_transfer_coeff = heat_transfer_coeff(external_state, external_nusselt_number, external_characteristic_length)
    
    return 1 / ((1 / (internal_heat_transfer_coeff * internal_area)) + (conductive_resistance) + (1 / (external_heat_transfer_coeff * external_area)))


def effectiveness_NTU_thermo(hot_stream_state, hot_stream_mass_flow, cold_stream_state, cold_stream_mass_flow, 
                             subtype, hot_inlet_temp, hot_outlet_temp, cold_inlet_temp, cold_outlet_temp, UA, n_shell_tube):
    hot_heat_capacity = hot_stream_state.Cp()
    cold_heat_capacity = cold_stream_state.Cp()
    return ht.effectiveness_NTU_method(mh=hot_stream_mass_flow, mc=cold_stream_mass_flow, 
                                Cph=hot_heat_capacity, Cpc=cold_heat_capacity,
                                subtype=subtype, 
                                Thi=hot_inlet_temp, Tho=hot_outlet_temp, 
                                Tci=cold_inlet_temp, Tco=cold_outlet_temp, 
                                UA=UA, n_shell_tube=n_shell_tube)
    
hot_state = flasher.flash(T=350, P=100000, zs=[0.7, 0.3])
cold_state = flasher.flash(T=300, P=100000, zs=[0.3, 0.7])


nusselt_1 = ht_internal_conv(
                state=cold_state, 
                state_mass_flow=100, 
                pipe_characteristic_length=100, 
                pipe_length=20, 
                absolute_pipe_roughness=0.0015
                )
print(nusselt_1)

nusselt_2 = ht_submerged_coil(
            ambient_state=hot_state, 
            heat_exchange_surface_temperature=340, 
            outer_characteristic_length=0.11, 
            horizontal=False
            )
print(nusselt_2)

overall_heat_transfer_coeff_1 = overall_heat_transfer_coeff(
        internal_state=cold_state,
        internal_nusselt_number=nusselt_1,
        internal_characteristic_length=0.1, 
        internal_area=1,
        external_nusselt_number=nusselt_2,
        external_state=hot_state, 
        external_characteristic_length=0.11, 
        external_area=1, 
        layers=[(0.012, 400, 10)]
        )
print(overall_heat_transfer_coeff_1)
#remember: [[thickness, thermal conductivity, area], ...]

#Without set UA:
print(effectiveness_NTU_thermo(
    hot_stream_state=hot_state, 
    hot_stream_mass_flow=100, 
    cold_stream_state=cold_state, 
    cold_stream_mass_flow=100, 
    subtype="counterflow", 
    hot_inlet_temp=hot_state.T, 
    hot_outlet_temp=None, 
    cold_inlet_temp=cold_state.T, 
    cold_outlet_temp=None, 
    UA=overall_heat_transfer_coeff_1, 
    n_shell_tube=None, 
))

#With set UA:
print(effectiveness_NTU_thermo(
    hot_stream_state=hot_state, 
    hot_stream_mass_flow=100, 
    cold_stream_state=cold_state, 
    cold_stream_mass_flow=100, 
    subtype="counterflow", 
    hot_inlet_temp=hot_state.T, 
    hot_outlet_temp=None, 
    cold_inlet_temp=cold_state.T, 
    cold_outlet_temp=None, 
    UA=10000, 
    n_shell_tube=None, 
))

def log_mean_temperature_difference(T_hot_in, T_hot_out, T_cold_in, T_cold_out):
    T1 = abs(T_hot_in - T_cold_in)
    T2 = abs(T_hot_out - T_cold_out)
    return (T1 - T2) / (np.log(T1 / T2))

def heat_duty(heat_transfer_coefficient, inner_characteristic_length, pipe_length, T_hot_in, T_hot_out, T_cold_in, T_cold_out):
    area = (math.pi * inner_characteristic_length) * pipe_length
    log_mean_temp_diff = log_mean_temperature_difference(T_hot_in, T_hot_out, T_cold_in, T_cold_out)
    return heat_transfer_coefficient * log_mean_temp_diff
