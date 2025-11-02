using Unitful
#using Symbolics, Nemo, Groebner
using Clapeyron
#using NLsolve

#using CoolProp : cp
import CoolProp as Cp



function solar_heat(sun_watts_per_m2::typeof(1.0u"W/m^2"), collector_area::typeof(1.0u"m^2"), collector_efficiency)
    return sun_watts_per_m2 * collector_area * collector_efficiency
end


function required_area_for_wattage(solar_heat::typeof(1.0u"W"), sun_watts_per_m2::typeof(1.0u"W/m^2"), collector_efficiency)
    return solar_heat / (sun_watts_per_m2 * collector_efficiency)
end

# Humidity Ratio (ω) is the ratio between the mass of water vapor in the air compared to the mass of the air
function humidity_ratio(partial_pressure_water::typeof(1.0u"Pa"))
    atmospheric_pressure = 101325u"Pa"
    #note that 0.622 is the ratio of the molar mass of water to the molar mass of dry air
    return uconvert(u"kg/kg", 0.622 * partial_pressure_water / (atmospheric_pressure - partial_pressure_water))
end

water_model = PR(["water"])

function relative_humidity(partial_pressure_water::typeof(1.0u"Pa"), air_temperature::typeof(1.0u"K"))
    water_p_sat_in_air = saturation_pressure(water_model, ustrip(air_temperature))
    return partial_pressure_water / water_p_sat_in_air
end

function calc_water_partial_pressure(relative_humidity, air_temperature::typeof(1.0u"K"))
    water_p_sat_in_air = saturation_pressure(water_model, air_temperature)[1]
    return relative_humidity * water_p_sat_in_air
end

#calc_water_partial_pressure(0.29, (25.0+273.13)u"K")
#saturation_pressure(water_model, (25.0+273.13)u"K")

kg_water_to_kg_air = humidity_ratio((calc_water_partial_pressure(0.9, 25.0u"°C" |> u"K")))

target_kg_s_of_water = uconvert(u"kg/s", (5u"kg" / 12u"hr"))

kg_s_of_process_air = target_kg_s_of_water / kg_water_to_kg_air 



# Something here for moist air enthalpy
# Define your models ONCE at the top of your script
air_components = ["nitrogen", "oxygen", "argon"] # Simplified dry air
water_and_air_components = ["water", "nitrogen", "oxygen", "argon"]
air_model = PR(air_components)
moist_air_model = PR(water_and_air_components)

function moist_air_enthalpy(T::typeof(1.0u"K"), P::typeof(1.0u"Pa"), humidity_ratio)
    # T and P are the temperature and total pressure of the moist air
    # ω is the humidity ratio (kg_w / kg_da)
    
    # Convert humidity ratio to mole fractions for Clapeyron
    water_molar_mass = 18.015e-3 # kg/mol
    dry_air_molar_mass = 28.96e-3 # kg/mol
    
    water_moles = humidity_ratio / water_molar_mass
    dry_air_moles = 1.0 / dry_air_molar_mass # basis of 1 kg dry air
    
    water_molar_concentration = moles_w / (moles_w + moles_da)
    air_molar_concentration = 1 - water_molar_concentration
    
    # Mole fractions of components in the total mixture
    # Assuming standard air composition (78% N2, 21% O2, 1% Ar)
    z = [water_molar_concentration, air_molar_concentration*0.78, air_molar_concentration*0.21, air_molar_concentration*0.01]
    
    # Calculate molar enthalpy (J/mol) of the mixture
    moist_air_molar_enthalpy = enthalpy(moist_air_model, ustrip(P), ustrip(T), z) # Gas phase
    
    # Convert to specific enthalpy (J/kg of dry air)
    total_moles = water_moles + dry_air_moles
    total_mass_mixture = humidity_ratio + 1.0 # kg
    h_specific_mixture = (moist_air_molar_enthalpy * total_moles) / total_mass_mixture # J/kg_mixture
    
    # This is tricky, often the enthalpy is desired per kg of DRY air.
    # For that we need enthalpy of dry air and enthalpy of water vapor separately.
    # Let's stick with the empirical one for now to keep it simple, but this is the rigorous way.
    # The empirical formula is fine for most applications near atmospheric pressure.
end


# Three below are just for mass balancing
function water_removal_mass_flow(dry_air_through_absorber_mass_flow::typeof(1.0u"kg/s"), humidity_ratio_in_absorber, humidity_ratio_out_absorber)
    return dry_air_through_absorber_mass_flow * (humidity_ratio_in_absorber - humidity_ratio_out_absorber)
end


function water_collector_mass_flow(dry_air_through_regenerator_mass_flow::typeof(1.0u"kg/s"), humidity_ratio_in_regenerator, humidity_ratio_out_regenerator)
    return dry_air_through_regenerator_mass_flow * (humidity_ratio_in_regenerator - humidity_ratio_out_regenerator)
end


# Strong desiccant enters absorber, weak desiccant leaves absorber
# Remember that by subtracting the mass flow of the weak desiccant that has absorbed water by the strong desiccant you get the mass flow of water removed
# weak_desiccant_mass_flow - strong_desiccant_mass_flow = water_removal_mass_flow
# The mass fraction of each determines how concentrated they are compared to water so strong would be like 0.9 and weak would be like 0.3
function desiccant_solution_concentration(strong_desiccant_mass_flow::typeof(1.0u"kg/s"), strong_desiccant_mass_fraction, 
                                            weak_desiccant_mass_flow::typeof(1.0u"kg/s"), weak_desiccant_mass_fraction)
    #NOTE: strong_desiccant_mass_flow * strong_desiccant_mass_fraction = weak_desiccant_mass_flow * weak_desiccant_mass_fraction
    mass_salt = m_s_strong * salt_mass_fractiontrong
    
    # The mass of the weak solution is the strong solution plus the water it absorbed
    m_s_weak = m_s_strong + m_w_removed
    
    # The concentration of the weak solution is the same mass of salt in a larger total mass
    weak_salt_mass_fraction = mass_salt / m_s_weak
    
    return (m_s_weak, weak_salt_mass_fraction)
end



# Note that strong equals concentrated desiccant and weak equals unconcentrated desiccant
function cooling_energy_balance(dry_air_through_absorber_mass_flow::typeof(1.0u"kg/s"), 
                                absorber_air_enthalpy_in::typeof(1.0u"J/kg"), absorber_air_enthalpy_out::typeof(1.0u"J/kg"), 
                                strong_desiccant_mass_flow::typeof(1.0u"kg/s"), strong_desiccant_enthalpy::typeof(1.0u"J/kg"), 
                                weak_desiccant_mass_flow::typeof(1.0u"kg/s"), weak_desiccant_enthalpy::typeof(1.0u"J/kg"))
    return dry_air_through_absorber_mass_flow * (absorber_air_enthalpy_in - absorber_air_enthalpy_out) + strong_desiccant_mass_flow * strong_desiccant_enthalpy - weak_desiccant_mass_flow * weak_desiccant_enthalpy
end

function heating_energy_balance(dry_air_through_regenerator_mass_flow::typeof(1.0u"kg/s"), 
                                reg_air_enthalpy_in::typeof(1.0u"J/kg"), reg_air_enthalpy_out::typeof(1.0u"J/kg"), 
                                strong_desiccant_mass_flow::typeof(1.0u"kg/s"), strong_desiccant_enthalpy::typeof(1.0u"J/kg"), 
                                weak_desiccant_mass_flow::typeof(1.0u"kg/s"), weak_desiccant_enthalpy::typeof(1.0u"J/kg"))
    return dry_air_through_regenerator_mass_flow * (reg_air_enthalpy_out - absorber_air_enthalpy_in) + strong_desiccant_mass_flow * strong_desiccant_enthalpy - weak_desiccant_mass_flow * weak_desiccant_enthalpy
end

water_model = PR(["water"])


function latent_heat_removed_from_air(water_removal_mass_flow::typeof(1.0u"kg/s"), air_temperature)
    water_model = ["water"]
    water_h_vap = enthalpy_vap(water_model, air_temperature)
    water_molar_mass = 18.015e-3u"kg/mol" 
    return water_removal_mass_flow * (water_h_vap * water_molar_mass)
end

function COP(latent_heat_removed_from_air::typeof(1.0u"kg/s"), solar_heat::typeof(1.0u"J/kg"))
    latent_heat_removed_from_air / solar_heat
end


mixture_model = PR(["water"])

T_solution = 300.0 # Temperature of the solution, e.g., 300 K (~27°C)
z_solution = [0.85, 0.15] # Composition: 85% water, 15% CaCl by mole fraction

# This tells you the total pressure above the liquid at equilibrium.
(P_bubble, V_l, V_v, y_vapor) = bubble_pressure(mixture_model, T_solution, z_solution)

# 4. Calculate the partial pressure of water
# The vapor composition is given by y_vapor = [y_water, y_CaCl]
# Since CaCl is a salt, its vapor pressure is negligible (y_CaCl ≈ 0).
# The total bubble pressure is essentially all due to water.
water_partial_pressure_sol = y_vapor[1] * P_bubble




# overall_mass_transfer_coefficient (empirical/to be found) describes how quickly water vapor can move from the air to the desiccant fluid
# specific_interfacial_packing (found from data sheet from a packing media provider) how much surface area a packing media provides per unit volume
function height_of_transfer_unit(air_mass_flow_rate::typeof(1.0u"kg/s"), column_cross_sectional_area::typeof(1.0u"m^2"), 
                                overall_mass_transfer_coefficient::typeof(1.0u"mol/(m^2*s)"), specific_interfacial_packing::typeof(1.0u"m^2/m^3"))
    air_molar_flow_rate = air_mass_flow_rate / 28.96u"g/mol"
    gas_molar_flux = air_molar_flow_rate / column_cross_sectional_area
    return gas_molar_flux / (overall_mass_transfer_coefficient * specific_interfacial_packing)
end

#enter humidity
start_relative_humidity = 0.90

end_relative_humidity = 0.50
#Can't go below 0.1 because that's lower than the desiccant's ceiling

kg_water_to_kg_air = humidity_ratio((calc_water_partial_pressure(start_relative_humidity, 25.0u"°C" |> u"K")))

target_kg_s_of_water = uconvert(u"kg/s", (5u"kg" / 12u"hr"))

kg_s_of_process_air = target_kg_s_of_water / kg_water_to_kg_air 

kg_s_of_process_air

required_HTU = height_of_transfer_unit(kg_s_of_process_air, 0.001u"m^2", 0.004u"mol/(m^2*s)", 200.0u"m^2/m^3")
# Search: typical overall mass transfer coefficient mol/(m^2*s) for calcium chloride desiccant with humid air
uconvert(u"m", required_HTU) 
#For above^ I got 0.22 m




function vapor_pressure_LiCl_solution(T_s::typeof(1.0u"K"), salt_mass_fraction)
    # --- Critical constants for water ---
    T_crit = 647.096u"K"
    P_crit = 22064.0u"kPa"

    # --- Constants for the polynomial fits of coefficients n1 to n6 ---
    # These are taken directly from the ASHRAE tables.
    A = [-7.859517, 1.84408, -11.7866, 22.6807, -15.9618, 1.801225]
    B = [-3.6666, -11.9332, 57.5822, -133.013, 114.333, -37.5953]
    C = [19.2666, 63.6332, -370.833, 854.166, -737.5, 241.666]
    D = [-25.0, -100.0, 500.0, -1166.66, 1000.0, -333.333]
    
    # --- Intermediate variables ---
    T = T_s # Use the input temperature
    θ = 1.0 - (T / T_crit)
    
    # --- Calculate coefficients n1 through n6 ---
    n = zeros(6)
    for i in 1:6
        n[i] = A[i] + B[i]*salt_mass_fraction + C[i]*salt_mass_fraction^2 + D[i]*salt_mass_fraction^3
    end
    
    # --- Calculate the summation term in the main equation ---
    sum_term = n[1]*θ + n[2]*θ^1.5 + n[3]*θ^3 + n[4]*θ^3.5 + n[5]*θ^4 + n[6]*θ^7.5
    
    # --- Calculate the final vapor pressure ---
    exponent = (T_crit / T) * sum_term
    P_v_sol = P_crit * exp(exponent)
    
    return uconvert(u"Pa", P_v_sol) # Ensure output is in Pascals
end

function calculate_y_eq_empirical(T_s::typeof(1.0u"K"), salt_mass_fraction, P_atm::typeof(1.0u"Pa"))
    P_v_sol = vapor_pressure_LiCl_solution(T_s, salt_mass_fraction)
    y_eq = P_v_sol / P_atm
    return y_eq
end


function calculate_NTU(model, y_in, y_out, salt_temperature_in, salt_temperature_out, salt_mass_fraction_in, salt_mass_fraction_out, atmospheric_pressure; num_steps=50)
    
    # Create the 'y' points for our slices, from outlet to inlet
    y_points = range(y_out, y_in, length=num_steps + 1)
    
    total_NTU = 0.0

    # Loop through each slice of the column
    for i in 1:num_steps
        # Get the air humidity at the start and end of this small slice
        y1 = y_points[i]
        y2 = y_points[i+1]
        
        # --- Find the corresponding desiccant properties at points y1 and y2 ---
        # We assume a linear relationship between the change in air humidity and the
        # change in desiccant properties. This is a standard and effective simplification.
        
        # Fraction of the way through the column (in terms of humidity change)
        frac1 = (y1 - y_out) / (y_in - y_out)
        frac2 = (y2 - y_out) / (y_in - y_out)

        # Linearly interpolate the desiccant temperature and concentration
        # Note the "cross-over": y_out (top) corresponds to salt_temperature_in (top)
        T_s1 = salt_temperature_out + (salt_temperature_in - salt_temperature_out) * frac1
        salt_mass_fraction1 = salt_mass_fraction_out + (salt_mass_fraction_in - salt_mass_fraction_out) * frac1

        T_s2 = salt_temperature_out + (salt_temperature_in - salt_temperature_out) * frac2
        salt_mass_fraction2 = salt_mass_fraction_out + (salt_mass_fraction_in - salt_mass_fraction_out) * frac2

        # --- Calculate the driving force at each end of the slice ---
        water_mole_fraction_at_eq1 = calculate_y_eq_empirical(T_s1, salt_mass_fraction1, atmospheric_pressure)
        water_mole_fraction_at_eq2 = calculate_y_eq_empirical(T_s2, salt_mass_fraction2, atmospheric_pressure)

        # This is the value we are integrating: 1 / (y - water_mole_fraction_at_eq)
        integrand1 = 1 / (y1 - water_mole_fraction_at_eq1)
        integrand2 = 1 / (y2 - water_mole_fraction_at_eq2)
        
        # --- Apply the Trapezoidal Rule for this one slice ---
        # dNTU = average_height * width
        avg_integrand = (integrand1 + integrand2) / 2.0
        dy = y2 - y1 # The "width" of our slice
        
        dNTU = avg_integrand * dy
        
        # Add the result for this slice to our running total
        total_NTU += dNTU
    end
    
    return total_NTU
end



# --- 2. Define Operating Conditions ---
atmospheric_pressure = 101325.0u"Pa"

# Air conditions (for a hot, humid day)
air_temperature_in = 24.0u"°C" |> u"K"
RH_in = start_relative_humidity # 90% Relative Humidity
# Desired outlet air conditions
air_temperature_out = 26.0u"°C" |> u"K" # Air heats up slightly
RH_out = end_relative_humidity # 40% Relative Humidity

# Desiccant conditions (counter-current)
# Strong, cool desiccant enters the top (where dry air exits)
salt_temperature_in = 22.0u"°C" |> u"K"
salt_mass_fraction_in = 0.45 # 45% CaCl2 by mass (strong)

# Weak, warmer desiccant leaves the bottom (where humid air enters)
salt_temperature_out = 27.0u"°C" |> u"K"# Heats up from absorbing water
salt_mass_fraction_out = 0.40 # 40% CaCl2 by mass (diluted)

# --- 3. Convert RH to Mole Fraction (y) ---
# We need a simple water model for this conversion
water_model = PR(["water"])
(p_sat_in, _, _) = saturation_pressure(water_model, ustrip(air_temperature_in))
(p_sat_out, _, _) = saturation_pressure(water_model, ustrip(air_temperature_out))

water_partial_pressure_in = RH_in * p_sat_in * u"Pa"
water_partial_pressure_out = RH_out * p_sat_out * u"Pa"

y_in = water_partial_pressure_in / atmospheric_pressure
y_out = water_partial_pressure_out / atmospheric_pressure

println("Air Inlet Mole Fraction (y_in): ", round(y_in, digits=4))
println("Air Outlet Mole Fraction (y_out): ", round(y_out, digits=4))
println("-"^30)

# --- 4. Run the Calculation ---
NTU_result = calculate_NTU(mixture_model, y_in, y_out, 
                            salt_temperature_in, salt_temperature_out, 
                            salt_mass_fraction_in, salt_mass_fraction_out, 
                            atmospheric_pressure; num_steps=50)
#For above^ I got 0.687 m
println("Calculated Number of Transfer Units (NTU): ", round(NTU_result, digits=2))

z = NTU_result * required_HTU

uconvert(u"m", z) 
#For above^ I got 0.15 m


function vapor_pressure_LiCl_solution(T_s::typeof(1.0u"K"), salt_mass_fraction)
    # --- Critical constants for water ---
    T_crit = 647.096u"K"
    P_crit = 22064.0u"kPa"

    # --- Constants for the polynomial fits of coefficients n1 to n6 ---
    # These are taken directly from the ASHRAE tables.
    A = [-7.859517, 1.84408, -11.7866, 22.6807, -15.9618, 1.801225]
    B = [-3.6666, -11.9332, 57.5822, -133.013, 114.333, -37.5953]
    C = [19.2666, 63.6332, -370.833, 854.166, -737.5, 241.666]
    D = [-25.0, -100.0, 500.0, -1166.66, 1000.0, -333.333]
    
    # --- Intermediate variables ---
    T = T_s # Use the input temperature
    θ = 1.0 - (T / T_crit)
    
    # --- Calculate coefficients n1 through n6 ---
    n = zeros(6)
    for i in 1:6
        n[i] = A[i] + B[i]*salt_mass_fraction + C[i]*salt_mass_fraction^2 + D[i]*salt_mass_fraction^3
    end
    
    # --- Calculate the summation term in the main equation ---
    sum_term = n[1]*θ + n[2]*θ^1.5 + n[3]*θ^3 + n[4]*θ^3.5 + n[5]*θ^4 + n[6]*θ^7.5
    
    # --- Calculate the final vapor pressure ---
    exponent = (T_crit / T) * sum_term
    P_v_sol = P_crit * exp(exponent)
    
    return uconvert(u"Pa", P_v_sol) # Ensure output is in Pascals
end

function calculate_y_eq_empirical(T_s::typeof(1.0u"K"), salt_mass_fraction, P_atm::typeof(1.0u"Pa"))
    P_v_sol = vapor_pressure_LiCl_solution(T_s, salt_mass_fraction)
    Cp.PropsSI("P", "T", T_s, "Q", 0, "INCOMP::MLI[0.21]")
    y_eq = P_v_sol / P_atm
    return y_eq
end

#Cp.PropsSI("P", "T", 300, "Q", 0, "INCOMP::MLI[0.1]")
#above returns 0 for some reason

#NOTE look into a direct solar falling film regenerator
#Just look into falling film absorbers and regenerators in general

#but not permeable to the liquid desiccant and if you put a black plate close enough so that the liquid is pulled towards it via surface tension, could the water boil out of the desiccant and go through the membrane without having the chance of reabsorbing in the desiccant

#in this case, y_eq_in is the mole fraction of water in the air that would be in equilibrium with the strong desiccant entering the absorber
#and y_eq_out is the mole fraction of water in the air that would be in equilibrium with the weak desiccant leaving the absorber
function log_mean_driving_force(y_in, y_out, y_eq_in, y_eq_out)
    atmospheric_pressure = 101325.0u"Pa"
    return ((y_in - y_eq_out) - (y_out - y_eq_in)) / log((y_in - y_eq_out) / (y_out - y_eq_in))
end

function overall_mass_transfer_coefficient(y_eq_in, y_eq_out, salt_mass_fraction_in, salt_mass_fraction_out, gas_side_mass_transfer_coefficient, liquid_side_mass_transfer_coefficient)
    delta_y_eq = y_eq_out - y_eq_in 
    
    moles_water_in = (1 - salt_mass_fraction_in) / 18.015
    moles_salt_in  = salt_mass_fraction_in / 42.39 #change later
    x_in = moles_water_in / (moles_water_in + moles_salt_in)
    moles_water_out = (1 - salt_mass_fraction_out) / 18.015
    moles_salt_out  = salt_mass_fraction_out / 42.39 #change later
    x_out = moles_water_out / (moles_water_out + moles_salt_out) 

    dx = abs(x_in - x_out)

    slope_of_equilbrium_line = delta_y_eq / dx
    #actual equation is 1 / overall_mass_transfer_coefficient = (1 / gas_side_mass_transfer_coefficient + slope_of_equilbrium_line / liquid_side_mass_transfer_coefficient)
    return 1 / (1 / gas_side_mass_transfer_coefficient + slope_of_equilbrium_line / liquid_side_mass_transfer_coefficient)
end

#Cp.PropsSI("D", "T", 300, "P", 100000, "INCOMP::MLI[0.2]")


function volumetric_flow_rate_from_overall_mass_transfer_coefficient_y_eq(overall_mass_transfer_coefficient, plate_area, log_mean_driving_force)
    return overall_mass_transfer_coefficient * plate_area * log_mean_driving_force
    #Note that y_eq can be subbed out for partial pressures instead in the log mean driving force
end

function gas_side_mass_transfer_coefficient(sherwood_number, diffusion_coefficient, characteristic_length)
    return (sherwood_number * diffusion_coefficient) / characteristic_length
    #sherwood_number = (gas_side_mass_transfer_coefficient * characteristic_length) / diffusion_coefficient
    #rearrange above for gas_side_mass_transfer_coefficient
    
    #diffusion_coefficient is for water vapor in air - probably 0.21cm^2/s
    
    #characteristic_length would be the length of the plate
end

function liquid_side_mass_transfer_coefficient(sherwood_number, diffusion_coefficient, film_thickness)
    return (sherwood_number * diffusion_coefficient) / film_thickness
    #sherwood_number = (liquid_side_mass_transfer_coefficient * film_thickness) / diffusion_coefficient
    #rearrange above for gas_side_mass_transfer_coefficient
    
    #diffusion_coefficient is for water vapor in air - probably 0.21cm^2/s

    #film thickness can be found from fluid properties and flow rate based on Nusselt's film theory for laminar flow 
end


#academic literature for "mass transfer correlations for falling film absorption" 
#or "Sherwood number for gas flow over a flat plate" and 
#"Sherwood number for falling liquid film."



function reynolds_number_with_coolprop(coolprop_fluid::String, temperature::typeof(1.0u"K"), pressure::typeof(1.0u"Pa"), 
                                        fluid_velocity::typeof(1.0u"m/s"), characteristic_length::typeof(1.0u"m"))
    density = Cp.PropsSI("D", "T", temperature, "P", pressure, coolprop_fluid)
    dynamic_viscosity = Cp.PropsSI("V", "T", temperature, "P", pressure, coolprop_fluid)
    return (density * fluid_velocity * characteristic_length) / (dynamic_viscosity)
end

function falling_film_reynolds_number(coolprop_fluid::String, temperature, pressure, fluid_velocity, plate_width, film_thickness)
    volumetric_flow = fluid_velocity * plate_width * film_thickness
    mass_flow = volumetric_flow * Cp.PropsSI("D", "T", temperature, "P", pressure, coolprop_fluid)
    mass_flow_per_unit_width = mass_flow / plate_width
    dynamic_viscosity = Cp.PropsSI("V", "T", temperature, "P", pressure, coolprop_fluid)
    return 4 * mass_flow_per_unit_width / dynamic_viscosity
end

function schmidt_number_with_coolprop(coolprop_fluid, temperature, pressure, mass_diffusivity)
    density = Cp.PropsSI("D", "T", temperature, "P", pressure, coolprop_fluid)
    dynamic_viscosity = Cp.PropsSI("V", "T", temperature, "P", pressure, coolprop_fluid)
    kinematic_viscosity = dynamic_viscosity / density
    return kinematic_viscosity / mass_diffusivity
end

function sherwood_flat_laminar_local(Re, Sc)
    #for 0.6 < Sc < 50
    return 0.332*Re^(1/2) * Sc^(1/3)
end

function sherwood_flat_laminar_average(Re, Sc)
    #for 0.6 < Sc < 50
    return 0.664*Re^(1/2) * Sc^(1/3)
end

function sherwood_turbulent_local(Re, Sc)
    #for 0.6 < Sc < 50
    #for Re_x < 10^8
    0.0296*Re^(4/5) * Sc^(1/3)
end

function sherwood_mixed_average(Re, Sc)
    #for 0.6 < Sc < 50
    #for Re_x < 10^8
    #for Re_xe < 5*10^5
    (0.037*Re^(4/5) - 871) * Sc^(1/3)
end

function sherwood_liquid(Re, Sc, film_thickness, characteristic_length)
    return (2 / sqrt(π)) * (Re * Sc * film_thickness / characteristic_length)^(1/2)
end
air_fluid = "air"

water_fluid = "INCOMP::MLI[0.24]"

temperature = 22.0u"°C" |> u"K"

pressure = 101325.0u"Pa"

fluid_velocity = 2.5u"m/s"

liquid_film_velocity = 0.3u"m/s"

characteristic_length = 0.5u"m"

plate_width = 1u"m"

plate_area = characteristic_length * plate_width 

water_in_salt_sol_mass_diffusivity = 1.6e-9u"m^2/s"

water_vapor_in_air_mass_diffusivity = 2.1e-5u"m^2/s"

film_thickness = 0.0005u"m"

atmospheric_pressure = 101325.0u"Pa"

volumetric_flow = liquid_film_velocity * plate_width * film_thickness
mass_flow = volumetric_flow * Cp.PropsSI("D", "T", temperature, "P", pressure, water_fluid)


gas_reynolds_number = reynolds_number_with_coolprop(
    air_fluid, 
    temperature, 
    pressure, 
    fluid_velocity, 
    characteristic_length
    )

gas_schmidt_number = schmidt_number_with_coolprop(
    air_fluid, 
    temperature, 
    pressure, 
    water_vapor_in_air_mass_diffusivity
    )

gas_sherwood_number = ustrip(sherwood_flat_laminar_average(ustrip(gas_reynolds_number), ustrip(gas_schmidt_number)))

gas_side_mass_transfer_coefficient_val = gas_side_mass_transfer_coefficient(
    gas_sherwood_number,
    water_vapor_in_air_mass_diffusivity,
    characteristic_length
)
#I got 0.007172847955176885 m s^-1

liquid_reynolds_number = falling_film_reynolds_number(
    water_fluid, 
    temperature, 
    pressure, 
    liquid_film_velocity, 
    plate_width, 
    film_thickness
    )

liquid_schmidt_number = schmidt_number_with_coolprop(
    water_fluid, 
    temperature, 
    pressure, 
    water_in_salt_sol_mass_diffusivity
    )

liquid_sherwood_number = ustrip(sherwood_liquid(ustrip(liquid_reynolds_number), ustrip(liquid_schmidt_number), film_thickness, characteristic_length))

liquid_side_mass_transfer_coefficient_val = liquid_side_mass_transfer_coefficient(
    liquid_sherwood_number,
    water_in_salt_sol_mass_diffusivity,
    film_thickness
)
#I got 0.0002018506017616128 m s^-1



RH_in = 0.9
RH_out = 0.45


air_temperature_in = 22.0u"°C" |> u"K"
air_temperature_out = 26.0u"°C" |> u"K"

salt_temperature_in = 22.0u"°C" |> u"K"
salt_temperature_out = 27.0u"°C" |> u"K"

salt_mass_fraction_in = 0.45
salt_mass_fraction_out = 0.40

water_model = PR(["water"])
(p_sat_in, _, _) = saturation_pressure(water_model, air_temperature_in)
(p_sat_out, _, _) = saturation_pressure(water_model, air_temperature_out)

water_partial_pressure_in = RH_in * p_sat_in 
water_partial_pressure_out = RH_out * p_sat_out 

y_in = water_partial_pressure_in / atmospheric_pressure #I got 0.019781635085671087
y_out = water_partial_pressure_out / atmospheric_pressure #I got 0.012689244541145886

y_eq_in = calculate_y_eq_empirical(salt_temperature_in, salt_mass_fraction_in, atmospheric_pressure) #I got 0.002879729537159945
y_eq_out = calculate_y_eq_empirical(salt_temperature_out, salt_mass_fraction_out, atmospheric_pressure) #I got 0.007573835104017883

overall_mass_transfer_coefficient_val = overall_mass_transfer_coefficient(y_eq_in, y_eq_out, salt_mass_fraction_in, salt_mass_fraction_out, gas_side_mass_transfer_coefficient_val, liquid_side_mass_transfer_coefficient_val)

log_mean_driving_force_val = log_mean_driving_force(y_in, y_out, y_eq_in, y_eq_out) #I got 0.010964979114838794

overall_mass_transfer_coefficient_val #I got 0.0013088095522355069 m s^-1

plate_area
log_mean_driving_force_val #I got 0.010964979114838794

volumetric_flow = volumetric_flow_rate_from_overall_mass_transfer_coefficient_y_eq(overall_mass_transfer_coefficient_val, plate_area, log_mean_driving_force_val) 
#I got 7.175534702781923e-6 m^3 s^-1

average_partial_pressure = ((y_in + y_out) / 2) * atmospheric_pressure

water_vapor_density = Cp.PropsSI("D", "T", temperature, "P", average_partial_pressure, "water")

mass_flow = volumetric_flow * water_vapor_density

uconvert(u"g/s", mass_flow)

kg_per_hour = uconvert(u"g/hr", mass_flow) 

kg_per_12_hours = uconvert(u"g/hr", mass_flow) * 12




#So currently I'm working on a solar liquid desiccant dehimidifier/water collector. I'm trying to size the abosorber and regenerator and I'm using this code to do so. Can you tell me if any of the values I used or got seem unreasonable? This is for a falling film absorber


