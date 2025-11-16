import numpy as np

import math 

import scipy
from scipy import optimize

import cantera as ct

from thermo import ChemicalConstantsPackage, CEOSGas, PRMIX, CEOSLiquid, FlashVL, equilibrium, EquilibriumStream, EquilibriumState, eos_mix, Flash


water = ct.Water()
water.TP = 300, 100000

methane = ct.Methane()
#ct.Solution("gri30.yaml")
methane.TP = 100, 100000

    
#print(short_cut_distliation_easy([water, methane], [0.5, 0.5], 0.99, 0.99))

                
            
            
        
    
"""outputs
    real output fractions
    update cantera fluid properties
    number of stages
    stage spacing
    10 values ranging from the tallest column with the least diameter to the shortest colum with the most diameter with their respective volumes
    vapor flowrate
    distillate flowrate
    bottoms flowrate
    reflux ration
    feed stage location
"""




def distillation_MESH(flasher, feed_stage, feed_temp, feed_pressure, feed_molar_flow, feed_mole_fractions, number_of_stages, stage_heat_exchange, 
                     product_withdrawal_locations, product_molar_flow_rates,
                     column_pressures, reboiler_pressure, condenser_pressure):
    #stage heat exchange is a list that contains the wattages for each column 
    Nc = len(feed_mole_fractions)
    
    #excluding reboiler and condenser
    N_stages = number_of_stages
    N_stages_total = number_of_stages + 1
    N_stages_model = number_of_stages + 2
    
    """
    j = 0: Reboiler
    j = 1 to N_stages: stages
    j = N_trays + 1: Condenser (N_total_stages - 1)
    """
    
    liquid_product_flows_calc = np.zeros(N_stages_model)
    vapor_product_flows_calc = np.zeros(N_stages_model)
    
    
    M_res = np.zeros((N_stages_model, Nc)) # N mass balances per component (stages 0 to N + 1)
    #E_res = np.zeros((N_stages_model, Nc)) # N+1 equilibrium relations (stages 0 to N) - though flash handles implicitly
    S_liquid_res = np.zeros(N_stages_model)
    S_vapor_res = np.zeros((N_stages_model)) # N+1 sum relations (L and V) - often implicit in flash
    H_res = np.zeros(N_stages_model) # N+1 heat balances (stages 0 to N)
    
    if column_pressures == None:
        column_pressures = []   
        for j in range(N_stages_model):
            column_pressures.append(condenser_pressure * (j / N_stages_total) + reboiler_pressure * (1 - j / N_stages_total))
            
    withdrawal_point_map = {loc: i for i, loc in enumerate(product_withdrawal_locations)}
    
    def mesh_solver(guess_vars):
        feed_state = flasher.flash(zs=feed_mole_fractions, T=feed_temp, P=feed_pressure)
        T_guess = guess_vars[0:(N_stages_model)]
        L_guess = guess_vars[(N_stages_model): (N_stages_model)*2]
        V_guess = guess_vars[(N_stages_model)*2: (N_stages_model)*3]
        x_1D = guess_vars[(N_stages_model)*3:len(guess_vars)]
        x = x_1D.reshape((N_stages_model, Nc))
        y = np.zeros((N_stages_model, Nc))
        
        HL_stages = np.zeros(N_stages_model)
        HV_stages = np.zeros(N_stages_model)
        #every time a product is withdrawn from a stage, this is incremented by one
        
        for j in range(0, N_stages_model):
            x_j_norm = x[j] / np.sum(x[j]) if np.sum(x[j]) > 1e-9 else x[j] # Avoid division by zero
            state = flasher.flash(zs=x_j_norm, T=T_guess[j], P=column_pressures[j])
        
            for i in range(Nc):
                    #Equilibrium
                    #if state.phase == "VL":
                        #y[j][i] = (state.Ks(state))[i] * x[j][i]
                    #else:
                    y[j][i] = (state.Psats()[i] / column_pressures[j]) * x[j][i]
            
            if j in withdrawal_point_map:
                withdrawal_index = withdrawal_point_map[j]
                total_prod_flow_at_j = product_molar_flow_rates[withdrawal_index]
                print("WITHDRAWAL", j)
                if state.phase == 'L':
                    liquid_product_flows_calc[j] = total_prod_flow_at_j
                    vapor_product_flows_calc[j] = 0.0
                elif state.phase == 'V':
                    liquid_product_flows_calc[j] = 0.0
                    vapor_product_flows_calc[j] = total_prod_flow_at_j
                elif state.phase == 'VL':
                    liquid_product_flows_calc[j] = total_prod_flow_at_j * (1 - state.VF)
                    vapor_product_flows_calc[j] = total_prod_flow_at_j * state.VF
                else: # Fallback for unexpected or critical phases (e.g., supercritical)
                    liquid_product_flows_calc[j] = total_prod_flow_at_j * 0.5 # Default to 50/50 split
                    vapor_product_flows_calc[j] = total_prod_flow_at_j * 0.5
            
            
            if j == 0: #for reboilerstat
                print("reboiler")
                for i in range(Nc):
                    #Mass
                    M_res[j][i] = L_guess[1] * x[j + 1][i] - V_guess[j] * y[j][i] - L_guess[j] * x[j][i]
                    if j in withdrawal_point_map:
                        M_res[j][i] -= liquid_product_flows_calc[j] * x[j][i]
                        M_res[j][i] -= vapor_product_flows_calc[j] * y[j][i]

            elif j == N_stages_total: #for condenser
                print("condenser")
                for i in range(Nc):
                    #Mass
                    M_res[j][i] = V_guess[j-1] * y[j-1][i] - L_guess[j] * x[j][i] - V_guess[j] * y[j][i]
                    if j in withdrawal_point_map:
                        M_res[j][i] -= liquid_product_flows_calc[j] * x[j][i]
                        M_res[j][i] -= vapor_product_flows_calc[j] * y[j][i]
                    
            elif j == feed_stage: #for feed stage
                print("feed stage")
                for i in range(Nc):
                    #Mass
                    M_res[j][i] = (L_guess[j + 1] * x[j + 1][i]) + (V_guess[j - 1] * y[j - 1][i]) + feed_molar_flow * feed_mole_fractions[i] - (L_guess[j] * x[j][i]) - (V_guess[j] * y[j][i])
                    if j in withdrawal_point_map:
                        M_res[j][i] -= liquid_product_flows_calc[j] * x[j][i]
                        M_res[j][i] -= vapor_product_flows_calc[j] * y[j][i]
                    
            else: #for other stages 
                print("other stages")
                for i in range(Nc):
                    #Mass
                    M_res[j][i] = (L_guess[j + 1] * x[j + 1][i]) + (V_guess[j - 1] * y[j - 1][i]) - (L_guess[j] * x[j][i]) - (V_guess[j] * y[j][i])
                    if j in withdrawal_point_map:
                        M_res[j][i] -= liquid_product_flows_calc[j] * x[j][i]
                        M_res[j][i] -= vapor_product_flows_calc[j] * y[j][i]
            
            
            S_liquid_res[j] = np.sum(x[j]) - 1
            S_vapor_res[j] = np.sum(y[j]) - 1
            
            print("    tray number", j)
            print("    tray temperature", T_guess[j])
            print("    tray pressure", column_pressures[j])
            print("    Mole fractions", state.zs)
            print("    M_res", M_res[j])
            print("    S_liquid_res", S_liquid_res[j])
            print("    S_vapor_res", S_vapor_res[j])
            print("    H_res", H_res[j])
            print("    liquid mole fractions", x[j])
            print("    liquid flow guess", L_guess[j])
            print("    vapor mole fracitons", y[j])
            print("    vapor flow guess", V_guess[j])
            print("    current phase state", state.phase)
            if state.phase == "VL":
                print("    K values with EOS", state.Ks(state))
                K_values_with_p_sat = []
                for i in range(Nc):
                    K_values_with_p_sat.append(state.Psats()[i] / column_pressures[j])
                print("    K values with p_sat / p", K_values_with_p_sat)
            else:
                K_values_with_p_sat = []
                for i in range(Nc):
                    K_values_with_p_sat.append(state.Psats()[i] / column_pressures[j])
                print("    K values with p_sat / p", K_values_with_p_sat)
            print("")
                
            flow_enthalpy_balances = 0
            """for i in range(len(x[j])):
                bottom_liquid = flasher.flash(zs=x[j - 1], T=T_guess[j - 1], P=column_pressures[j - 1])
                top_vapor = flasher.flash(zs=x[j - 1], T=T_guess[j + 1], P=column_pressures[j + 1])
                current_state = flasher.flash(zs=x[j - 1], T=T_guess[j], P=column_pressures[j])
                #find the current column energy balance with the enthalpy of the bottom column liquid flowrate, top column vapor flowrate, current column liquid flowrate, 
                #and current column vapor flowrate * their respective enthalpies
                flow_enthalpy_balances += L_guess[j - 1] * bottom_liquid.LH + V_guess[j + 1] * top_vapor.VH - L_guess[j] * current_state.LH - V_guess[j] * current_state.VH 
                """
                
        for j in range(0, N_stages_model):
            liq_state = flasher.flash(zs=x[j], T=T_guess[j], P=column_pressures[j])
            print("vap state test", y[j], T_guess[j], column_pressures[j])
            vap_state = flasher.flash(zs=y[j], T=T_guess[j], P=column_pressures[j])
            try:
                HL_stages[j] = liq_state.liquid0.H()
                print("                              HL Success!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", liq_state.phase)
            except:
                HL_stages[j] = 0
                H_res[j] += 100000
                for i in range(Nc):
                    M_res[j][i] += 10000
                S_liquid_res[j] += 10000
                S_vapor_res[j] += 10000
                print("                              ZONKERS! HL Failiure!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", liq_state.phase)
            try:
                HV_stages[j] = vap_state.gas.H()
                print("                              HV Success!", "stage:", j, "y_j:", y[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", vap_state.phase)
            except:
                HV_stages[j] = 0
                H_res[j] += 100000
                for i in range(Nc):
                    M_res[j][i] += 10000
                S_liquid_res[j] += 10000
                S_vapor_res[j] += 10000
                print("                              ZOINKS! HV Failiure!", "stage:", j, "y_j:", y[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", vap_state.phase)
            """
            x_j_norm = x[j] / np.sum(x[j]) if np.sum(x[j]) > 1e-9 else x[j] # Avoid division by zero
            liq_state = flasher.flash(zs=x_j_norm, T=T_guess[j], P=column_pressures[j])
            y_j_norm = y[j] / np.sum(y[j]) if np.sum(y[j]) > 1e-9 else y[j] # Avoid division by zero
            vap_state = flasher.flash(zs=y_j_norm, T=T_guess[j], P=column_pressures[j])
            """
            
        for j in range(0, N_stages_model):
            #different H residuals for reboiler, condenser, and all other stages
            if j == 0: #Reboiler
                H_res[j] += L_guess[j + 1]*HL_stages[j + 1] - L_guess[j]*HL_stages[j] - V_guess[j]*HV_stages[j] + stage_heat_exchange[j]
                #H_res[j] = L_guess[j + 1]*HL_stages[j + 1] + V_guess[j - 1]*HV_stages[j - 1] - L_guess[j]*HL_stages[j] - V_guess[j]*HV_stages[j] + stage_heat_exchange[j]
            elif j == N_stages_total: #Condenser
                H_res[j] += V_guess[j-1]*HV_stages[j-1] - L_guess[j]*HL_stages[j] - V_guess[j]*HV_stages[j] + stage_heat_exchange[j]
                #H_res[j] = L_guess[j + 1]*HL_stages[j + 1] + V_guess[j - 1]*HV_stages[j - 1] - L_guess[j]*HL_stages[j] - V_guess[j]*HV_stages[j] + stage_heat_exchange[j]
            else: #Other Stages
                H_res[j] += L_guess[j + 1]*HL_stages[j + 1] + V_guess[j - 1]*HV_stages[j - 1] - L_guess[j]*HL_stages[j] - V_guess[j]*HV_stages[j] + stage_heat_exchange[j]
            
            #account for feed_stage enthalpy change if that stage is one of the feed_stages
            if j == feed_stage:
                H_res[j] += feed_molar_flow * feed_state.H()
            #accounting for enthalpy taken out when components are withdrawn as products
            if j in withdrawal_point_map:
                H_res[j] -= liquid_product_flows_calc[j] * HL_stages[j]
                H_res[j] -= vapor_product_flows_calc[j] * HV_stages[j]
        
            #now take the previous calculated energy balance and use the added enthalpy from the feed flow (if any) and the heat exchange for the stage (either losing heat or being heated up)
            
        return np.concatenate((M_res.flatten(), S_liquid_res.flatten(), S_vapor_res.flatten(), H_res.flatten()))
    
    # 1. Temperature Guesses (T_guess)
    T_guess = np.zeros(N_stages_model)
    # Estimate typical operating temperatures or use feed_temp as a guide
    # Reboiler is hotter, condenser is colder
    T_reboiler_initial = feed_temp + 20 # Example: 20K hotter than feed
    T_condenser_initial = feed_temp - 20 # Example: 20K colder than feed
    
    # Create a linear temperature profile from reboiler (j=0) to condenser (j=N_stages_total)
    T_guess = np.linspace(T_reboiler_initial, T_condenser_initial, N_stages_model)
    T_guess[T_guess < 1e-6] = 1e-6 # Ensure positive temperatures

    # 2. Flow Rate Guesses (L_guess, V_guess)
    L_guess = np.zeros(N_stages_model)
    V_guess = np.zeros(N_stages_model)

    F = feed_molar_flow # Total feed molar flow
    
    # Crude estimation of distillate (D) and bottoms (B) products based on feed flow
    # This is a very rough guess, adjust based on expected separation
    D_guess = F * 0.5 # Example: 50% of feed goes to distillate
    B_guess = F * 0.5 # Example: 50% of feed goes to bottoms

    # Assume a typical reflux ratio (R) and boilup ratio
    R_guess = 1.5 # Reflux Ratio (L/D)
    
    L_reflux_guess = R_guess * D_guess # Liquid reflux from condenser
    V_boilup_guess = (R_guess + 1) * D_guess # Vapor flow from reboiler (assuming V_boilup ~ (R+1)D)

    # Set boundary conditions for L and V guesses
    # Reboiler (j=0)
    L_guess[0] = B_guess # Liquid bottoms product
    V_guess[0] = V_boilup_guess # Vapor leaving reboiler

    # Condenser (j=N_stages_total)
    L_guess[N_stages_total] = L_reflux_guess # Liquid reflux to column
    V_guess[N_stages_total] = 1e-6 # Very small vapor flow out (assuming total condenser / very little non-condensables)

    # Linear interpolation for intermediate stages
    # L decreases from bottom to top, V increases from bottom to top
    L_guess_intermediate = np.linspace(L_guess[0], L_guess[N_stages_total], N_stages_model)
    V_guess_intermediate = np.linspace(V_guess[0], V_guess[N_stages_total], N_stages_model)
    
    # Adjust for product withdrawals if they significantly alter flows
    # For now, stick with simple linear profile, harder to account for side streams
    # without knowing their magnitude relative to total column flow.
    L_guess = L_guess_intermediate
    V_guess = V_guess_intermediate
    
    # Ensure flows are positive
    L_guess[L_guess < 1e-6] = 1e-6
    V_guess[V_guess < 1e-6] = 1e-6

    # 3. Liquid Mole Fractions Guesses (x)
    x = np.zeros((N_stages_model, Nc))

    # Assume components are ordered by volatility (e.g., component 0 is lightest, Nc-1 is heaviest)
    # This is a critical assumption for setting gradients.
    
    # Create a linear gradient for each component's mole fraction across the column
    for i in range(Nc):
        light_comp_bottom_frac = 0.05 # e.g., 5% of light component in bottoms
        light_comp_top_frac = 0.95    # e.g., 95% of light component in distillate

        heavy_comp_bottom_frac = 0.95 # e.g., 95% of heavy component in bottoms
        heavy_comp_top_frac = 0.05    # e.g., 5% of heavy component in distillate

        # Create end compositions for linear interpolation
        x_bottom_target = np.zeros(Nc)
        x_top_target = np.zeros(Nc)

        if Nc == 1:
            x_bottom_target[0] = 1.0
            x_top_target[0] = 1
        else:
            # Assign heaviest and lightest ends
            x_bottom_target[0] = light_comp_bottom_frac # Lightest comp in bottom
            x_bottom_target[Nc-1] = heavy_comp_bottom_frac # Heaviest comp in bottom
            x_top_target[0] = light_comp_top_frac # Lightest comp in top
            x_top_target[Nc-1] = heavy_comp_top_frac # Heaviest comp in top

            # Distribute remaining feed components proportionally to maintain sum=1
            remaining_sum_bottom = np.sum(feed_mole_fractions[1:-1])
            remaining_sum_top = np.sum(feed_mole_fractions[1:-1])

            if Nc > 2 and remaining_sum_bottom > 1e-9 and remaining_sum_top > 1e-9:
                x_bottom_target[1:-1] = feed_mole_fractions[1:-1] * (1.0 - x_bottom_target[0] - x_bottom_target[Nc-1]) / remaining_sum_bottom
                x_top_target[1:-1] = feed_mole_fractions[1:-1] * (1.0 - x_top_target[0] - x_top_target[Nc-1]) / remaining_sum_top
            else: # If only 2 components or feed is mostly end components
                 if Nc > 1: # If only 2 components
                    x_bottom_target[0] = 1.0 - x_bottom_target[Nc-1]
                    x_top_target[Nc-1] = 1.0 - x_top_target[0]
                 
        # Normalize the target compositions to ensure they sum to 1.0
        x_bottom_target /= np.sum(x_bottom_target)
        x_top_target /= np.sum(x_top_target)

        # Ensure no negative or extremely small values
        x_bottom_target[x_bottom_target < 1e-9] = 1e-9
        x_top_target[x_top_target < 1e-9] = 1e-9
        x_bottom_target /= np.sum(x_bottom_target) # Re-normalize after min check
        x_top_target /= np.sum(x_top_target) # Re-normalize after min check


        # Now, populate x with linear interpolation
        for j in range(N_stages_model):
            for i in range(Nc):
                # j=0 is reboiler (bottom), j=N_stages_total is condenser (top)
                x[j, i] = np.interp(j, [0, N_stages_total], [x_bottom_target[i], x_top_target[i]])
            
            # Ensure each stage's mole fractions sum to 1.0 after linear interpolation
            if np.sum(x[j]) > 1e-9:
                x[j, :] /= np.sum(x[j])
            else: # Fallback if sum is zero, make it even split (highly unlikely with sensible targets)
                x[j, :] = np.ones(Nc) / Nc

            # Ensure mole fractions are within [0, 1] range after normalization (due to interp)
            x[j][x[j] < 0] = 0 # Should not happen with np.interp if targets are 0-1
            x[j][x[j] > 1] = 1 # Should not happen
    
    guess_vars = np.concatenate((T_guess, L_guess, V_guess, x.flatten()))
    
    # --- Define bounds for scipy.optimize.least_squares ---
    num_vars = len(guess_vars)

    lb = np.zeros(num_vars)
    ub = np.full(num_vars, np.inf)

    current_idx = 0
    # Temperature bounds
    lb[current_idx : current_idx + N_stages_model] = 330
    ub[current_idx : current_idx + N_stages_model] = 400
    current_idx += N_stages_model

    # Liquid flow bounds
    lb[current_idx : current_idx + N_stages_model] = 1e-6
    ub[current_idx : current_idx + N_stages_model] = 1000 # Max reasonable flow rate
    current_idx += N_stages_model

    # Vapor flow bounds
    lb[current_idx : current_idx + N_stages_model] = 1e-6
    ub[current_idx : current_idx + N_stages_model] = 1000 # Max reasonable flow rate
    current_idx += N_stages_model

    
    lb[current_idx : current_idx + N_stages_model * (Nc)] = 1e-9
    ub[current_idx : current_idx + N_stages_model * (Nc)] = 1
    
    # Use least_squares with bounds for more stable optimization
    # Increased max_nfev, tightened tolerances for a better solution
    
    result = scipy.optimize.least_squares(mesh_solver, guess_vars, method="dogbox", bounds=(lb, ub), ftol=1e-6, xtol=1e-6, gtol=1e-6, max_nfev=20000, verbose=1)
    print("M_res")
    print(M_res)
    print(sum(M_res.flatten()))
    print("")
    
    print("S_liquid_res")
    print(S_liquid_res)
    print(sum(S_liquid_res.flatten()))
    print("")
    
    print("S_vapor_res")
    print(S_vapor_res)
    print(sum(S_vapor_res.flatten()))
    print("")
    
    print("H_res")
    print(H_res)
    print(sum(H_res.flatten()))
    print("")
    print("")
    
    print("T_guess")
    print(T_guess)
    print(sum(T_guess))
    print("")
    
    print("L_guess")
    print(L_guess)
    print(sum(L_guess))
    print("")
    
    print("V_guess")
    print(V_guess)
    print(sum(L_guess))
    print("")
    
    print("x_guess")
    print(x)
    print(sum(x.flatten()))
    print("")
    print("")
    print(result)

    # Combine all parts into guess_vars
    
    
    #return scipy.optimize.root(mesh_solver, guess_vars, method='lm', tol=1e-6)
    


constants, correlations = ChemicalConstantsPackage.from_IDs(['Acetic Acid', 'Water'])
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liq = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashVL(constants, correlations, liquid=liq, gas=gas)
mole_fractions = [0.5, 0.5]
test_flash = flasher.flash(zs=[0.5, 0.5], T=300, P=100000)
print("TEST", test_flash.T)
print("TEST", test_flash.Ks(test_flash))



feed_stage = 2
feed_temp = 273.13 + 100  # K (100°C, near benzene boiling point)
feed_pressure = 101325  # Pa (1 atm)
feed_molar_flow = 100  # mol/s

feed_mole_fractions = [0.5, 0.5]  
number_of_stages = 20

stage_adiabatic_specs = [0.0 for _ in range(number_of_stages)]
print(stage_adiabatic_specs)
stage_heat_exchange = ([530_000] + stage_adiabatic_specs + [-3_000_000])  # adiabatic stages
print(stage_heat_exchange)

product_withdrawal_locations = [0, 3]

product_molar_flow_rates=[45, 45]
    
reboiler_pressure = 60000  # Pa
condenser_pressure = 20000  # slight pressure drop

print(
    distillation_MESH(flasher=flasher, feed_stage=feed_stage, feed_temp=feed_temp, feed_pressure=feed_pressure, feed_molar_flow=feed_molar_flow, feed_mole_fractions=feed_mole_fractions, 
                 number_of_stages=number_of_stages, stage_heat_exchange=stage_heat_exchange, 
                 product_withdrawal_locations=product_withdrawal_locations, product_molar_flow_rates=product_molar_flow_rates,
                 column_pressures=None, 
                 reboiler_pressure=reboiler_pressure, condenser_pressure=condenser_pressure)
)