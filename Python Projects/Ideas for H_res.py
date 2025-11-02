"""
            liquid_stream_calculator = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
            liquid_stream_calculator.zs = x[j] # Set the composition
            liquid_stream_calculator.T = T_guess[j] # Set the temperature
            liquid_stream_calculator.P = column_pressures[j] # Set the pressure
            HL_stages[j] = liquid_stream_calculator.H() # Get the molar enthalpy
            print("                              HL Success!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", liquid_stream_calculator.phase)

            # Calculate HV_stages[j]: Enthalpy of the vapor stream (composition y[j])
            # Create a temporary gas-phase mixture object for the current stage's vapor composition
            # This assumes y[j] is the composition of the vapor stream flowing between stages
            vapor_stream_calculator = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
            vapor_stream_calculator.zs = y[j] # Set the composition
            vapor_stream_calculator.T = T_guess[j] # Set the temperature
            vapor_stream_calculator.P = column_pressures[j] # Set the pressure
            HV_stages[j] = vapor_stream_calculator.H() # Get the molar enthalpy
            print("                              HV Success!", "stage:", j, "y_j:", y[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", vapor_stream_calculator.phase)
            """
            ##
            #OR
            ##
    
            # The state object now contains the saturated liquid properties.
            # The check for 'if state_L.liquid' is robust in case the solver
            # goes to an invalid composition (e.g., negative fractions).
            if state_L.liquid0:
                HL_stages[j] = state_L.liquid0.H()
                print("                              HL Success!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "T_bubble", state_L.T, "Column Pressure:", column_pressures[j], "state:", state_L.phase)
            else:
                # This branch indicates a serious problem with the solver's state.
                # Handle it by setting a penalty or a safe value.
                HL_stages[j] = 1e6 # Or some other penalty value
                print("                              ZONKERS! HL Failiure!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "T_bubble", state_L.T, "Column Pressure:", column_pressures[j], "state:", future_liq_states[j].phase)
                print(f"CRITICAL WARNING: Could not find bubble point for stage {j}")

            # --- Calculate Saturated Vapor Enthalpy and Dew Temperature ---
            # We do the same for the vapor stream: find the dew point state (VF=1).
            
            if state_V.gas:
                HV_stages[j] = state_V.gas.H()
                print("                              HV Success!", "stage:", j, "v_j:", y[j], "T_guess:", T_guess[j], "T_bubble", state_V.T, "Column Pressure:", column_pressures[j], "state:", state_V.phase)
            else:
                HV_stages[j] = 1e6 # Or some other penalty value
                print(f"CRITICAL WARNING: Could not find dew point for stage {j}")
            
        """
            liq_state = flasher.flash(zs=x[j], T=T_guess[j], VF=0)
            vap_state = flasher.flash(zs=y[j], T=T_guess[j], VF=1)
            if liq_state.phase == "VL" or liq_state.phase == "L":
                HL_stages[j] = liq_state.liquid0.H()
                print("                              HL Success!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", liq_state.phase)
            else:
                HL_stages[j] = 0
                print("                              ZONKERS! HL Failiure!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", future_liq_states[j].phase)
            
            if vap_state.phase == "VL" or vap_state.phase == "V":
                HL_stages[j] = vap_state.gas.H()
                print("                              HL Success!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", vap_state.phase)
            else:
                HV_stages[j] = 0
                print("                              ZOINKS! HV Failiure!", "stage:", j, "y_j:", y[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", vap_state.phase)
        """
        """
            vap_state = flasher.flash(zs=y[j], T=T_guess[j], P=column_pressures[j])
            total_succeses = 0
            try:
                HL_stages[j] = future_liq_states[j].liquid0.H()
                total_succeses += 1
                print("                              HL Success!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", future_liq_states[j].phase)
            except:
                HL_stages[j] = 100000000000
                #OR
                #future_liq_states[j].H() * (1 - future_liq_states[j].quality)
                #OR
                #HL_stages[j] = 0
                print("                              ZONKERS! HL Failiure!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", future_liq_states[j].phase)
            try:
                HV_stages[j] = vap_state.gas.H()
                total_succeses += 1
                print("                              HV Success!", "stage:", j, "y_j:", y[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", vap_state.phase)
            except:
                HV_stages[j] = 100000000000
                #OR
                #vap_state.H() * vap_state.quality
                
                #HV_stages[j] = 0
                print("                              ZOINKS! HV Failiure!", "stage:", j, "y_j:", y[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", vap_state.phase)
        """
        #Note to future self: ask why the liquid state is always in the vapor state and the vapor state is always in the liquid state
        #OH WAIT! I think it's going towards this because we're giving it a higher enthalpy value whenever it does actually find a state in two phase
        #Maybe we should give the solver another set of variables (MESH + eval) that are either 1 or 0 depending on whether or not it reached a two phase state
        #the only issue is that if we append a huge value, it just ignores it and solves for another variable
        """
        ZONKERS! HL Failiure! stage: 0 x_j: [2.65933514e-04 9.32176503e-01] T_guess: 400.0762169289959 Column Pressure: 150000.0 state: V
        ZOINKS! HV Failiure! stage: 0 y_j: [9.30493858e-04 1.53088403e+00] T_guess: 400.0762169289959 Column Pressure: 150000.0 state: L
        ZONKERS! HL Failiure! stage: 1 x_j: [0.25541848 0.67618967] T_guess: 387.7919953570612 Column Pressure: 120000.00000000001 state: V
        ZOINKS! HV Failiure! stage: 1 y_j: [0.77106296 0.94225132] T_guess: 387.7919953570612 Column Pressure: 120000.00000000001 state: L
        ZONKERS! HL Failiure! stage: 2 x_j: [0.54102984 0.39152247] T_guess: 374.45570887595346 Column Pressure: 90000.0 state: V
        ZOINKS! HV Failiure! stage: 2 y_j: [1.40911104 0.46215256] T_guess: 374.45570887595346 Column Pressure: 90000.0 state: L
        ZONKERS! HL Failiure! stage: 3 x_j: [0.85465171 0.06797986] T_guess: 360.51147763857205 Column Pressure: 60000.0 state: V
        ZOINKS! HV Failiure! stage: 3 y_j: [2.0330087  0.07187323] T_guess: 360.51147763857205 Column Pressure: 60000.0 state: L
        """
        ##
        #OR
        ##
        """
        for j in range(0, N_stages_model):
            liq_state = flasher.flash(zs=x[j], T=T_guess[j], P=column_pressures[j])
            vap_state = flasher.flash(zs=y[j], T=T_guess[j], P=column_pressures[j])
            try:
                HL_stages[j] = liq_state.liquid0.H()
                print("                              HL Success!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", liq_state.phase)
            except:
                print(x[j], column_pressures[j])
                bubble_point_state = flasher.flash(zs=x[j], T=T_guess[j], VF=0.0001)
                HL_stages[j] = bubble_point_state.liquid0.H()
                print("use this", bubble_point_state.T)
                print("                              ZONKERS! HL Failiure!", "stage:", j, "x_j:", x[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", liq_state.phase)
            try:
                HV_stages[j] = vap_state.gas.H()
                print("                              HV Success!", "stage:", j, "y_j:", y[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", vap_state.phase)
            except:
                print(y[j], column_pressures[j])
                dew_point_state = flasher.flash(zs=x[j], T=T_guess[j], VF=0.9999)
                HV_stages[j] = dew_point_state.gas.H()
                print("use this", dew_point_state.T)
                print("                              ZOINKS! HV Failiure!", "stage:", j, "y_j:", y[j], "T_guess:", T_guess[j], "Column Pressure:", column_pressures[j], "state:", vap_state.phase)
        """