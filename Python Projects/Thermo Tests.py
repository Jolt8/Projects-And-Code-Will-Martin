from thermo import ChemicalConstantsPackage, CEOSGas, PRMIX, CEOSLiquid, FlashVL, equilibrium, EquilibriumStream, EquilibriumState
#constants = ChemicalConstantsPackage(names=['carbon dioxide', 'hexane'], CASs=['124-38-9', '110-54-3'], MWs=[44.0095, 86.17536], omegas=[0.2252, 0.2975], Pcs=[7376460.0, 3025000.0], Tbs=[194.67, 341.87], Tcs=[304.2, 507.6], Tms=[216.65, 178.075])
#correlations = PropertyCorrelationsPackage(constants=constants, skip_missing=True,
                                           #HeatCapacityGases=[HeatCapacityGas(poly_fit=(50.0, 1000.0, [-3.1115474168865828e-21, 1.39156078498805e-17, -2.5430881416264243e-14, 2.4175307893014295e-11, -1.2437314771044867e-08, 3.1251954264658904e-06, -0.00021220221928610925, 0.000884685506352987, 29.266811602924644])),
                                                              #HeatCapacityGas(poly_fit=(200.0, 1000.0, [1.3740654453881647e-21, -8.344496203280677e-18, 2.2354782954548568e-14, -3.4659555330048226e-11, 3.410703030634579e-08, -2.1693611029230923e-05, 0.008373280796376588, -1.356180511425385, 175.67091124888998]))])
constants, correlations = ChemicalConstantsPackage.from_IDs(['methane', 'water'])
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)
state = flasher.flash(zs=[0.5, 0.5], T=300, P=100000)
stream = EquilibriumStream(flasher=flasher, zs=[0.5, 0.5], T=300, P=100000, m=1)

test = flasher.flash(zs=[0.5, 0.5], H_mass=300, P=100000)
print(test.T)

test2 = flasher.flash(zs=[0.5, 0.5], H_mass=6000, P=100000)
print(test2.T)
print(test.Hvaps())

test2 = flasher.flash(zs=[0.5, 0.5], P=100000, VF=0.5)
print("looking for", test2.H())


def recalc(Equilibrium_Stream_obj:object):
    #store all values
    flasher=Equilibrium_Stream_obj.flasher
    zs=Equilibrium_Stream_obj.zs
    T=Equilibrium_Stream_obj.T
    P=Equilibrium_Stream_obj.P
    m=Equilibrium_Stream_obj.m
    
    return EquilibriumStream(flasher=flasher, zs=zs, T=T, P=P, m=m)
    #Equilibrium_Stream_obj.__init__(flasher=flasher, zs=zs, T=T, P=P, m=m, existing_flash=flasher) #= EquilibriumStream(flasher=flasher, zs=zs, T=T, P=P, m=m, existing_flash=None)

print(state.quality, state.T, state.P)

stream.T = 250 
stream.P = 100000
stream = recalc(stream) #idealy, this would happen every time an attribute was changed 
print(state.quality, state.T, state.P)

print(state.Ks(state, state.liquid0))
print(state.Ks(state, state.gas))
print(state.Ks(state, state.liquid0))

print(state.liquid0.H())
print(state.gas.H())
print("TETLJDLKFJ", state.liquid_bulk, state.LF, state.VF)

state1 = flasher.flash(zs=[0.00001, 1], T=500, P=100000)
try:
    print("state 1 liquid", state1.Ks(state1.liquid0))
except:
    print("state 1 gas", state1.Ks(state1.gas))
finally:
    print("state 1 all", state1.Ks(state1))


state2 = flasher.flash(zs=[1, 0.00001], T=500, P=100000)
try:
    print("state 2 liquid", state2.Ks(state2.liquid0))
except:
    print("state 2 gas", state2.Ks(state2.gas))
finally:
    print("state 2 all", state2.Ks(state2))
    
    
constants, correlations = ChemicalConstantsPackage.from_IDs(['methane', 'toluene', "benzene"])
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)
state = flasher.flash(zs=[0.5, 0.5, 0.5], T=300, P=100000)

state_for_ks = flasher.flash(zs=[0.356676, 0.27946676, 0.9], T=300.4633267057478, P=813333.33333333333)

print("state for Ks", state_for_ks.phase)
print("Ks calc all", state_for_ks.Ks(state_for_ks))
#print(dir(state_for_ks))
print("Ks calc", state_for_ks.Ks(state_for_ks.liquid0))







"""
def recalc(self):
    #store all values
    flasher=self.flasher
    zs=self.zs
    T=self.T
    P=self.P
    m=self.m
    
    #call it again 
    self = self.__init__(flasher=flasher, zs=zs, T=T, P=P, m=1)
    #self = EquilibriumStream(flasher=flasher, zs=zs, T=T, P=P, m=1)

EquilibriumStream.recalc = recalc
""" 
"""
state = flasher.flash(P=100000, T=310, zs=[0.5, 0.5])
print(state.quality)
state.T = 400
state = flasher.flash(P=state.P, T=state.T, zs=state.zs)
print(state.quality)
"""
"""
print("phase count", state.phase_count)
print("bulk cp", state.bulk.Cp())
print("flash specs", state.flash_specs)
print("melting temperature", state.Tms)
print("liquid enthalpy", state.liquid0.H())
print("gas enthalpy", state.gas.H())
print("quality", state.quality)
"""
