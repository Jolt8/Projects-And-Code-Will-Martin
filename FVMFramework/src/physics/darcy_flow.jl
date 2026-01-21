
function get_darcy_mass_flux(rho_avg, permeability, viscosity, pressure_a, pressure_b, area, dist) #not sure what k_average would represent in respect to diffusion
    p_grad = (pressure_b - pressure_a) / dist
    m_dot = -rho_avg * (permeability / viscosity) * p_grad * area
    return m_dot
end

function continuity_and_momentum_darcy(
    #NOTE!!:
    #this also returns face_m_dot even though it's a f!() function

    #cached vars
    mw_avg_cache_a, mw_avg_cache_b,
    #mutated vars
    du_pressure_a, du_pressure_b,
    #u values
    temp_a, temp_b,
    pressure_a, pressure_b,
    #geometry data
    area, norm, dist,
    vol_a, vol_b,
    #other props 
    rho_a, rho_b, #kinda a u value because it changes with time but not explicitly tracked through u values
    mu_a, mu_b,
    permeability,
    )

    rho_avg = 0.5 * (rho_a + rho_b)
    mu_avg = 0.5 * (mu_a + mu_b)

    face_m_dot = pressure_numerical_flux(rho_avg, permeability, mu_avg, pressure_a, pressure_b, area, dist)

    term_a = (face_m_dot / vol_a) * (R_gas * temp_a / mw_avg_cache_a)
    du_pressure_a -= term_a

    term_b = (face_m_dot / vol_b) * (R_gas * temp_b / mw_avg_cache_b)
    du_pressure_b += term_b

    return face_m_dot
end

