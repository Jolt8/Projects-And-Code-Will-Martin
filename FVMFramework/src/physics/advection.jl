#use get_cell_rho from helper functions to get rho for cell a and b

function species_advection!(
    #mutated vars
    du_species_mass_fractions_a, du_species_mass_fractions_b,
    #u values
    face_m_dot,
    species_mass_fractions_a, species_mass_fractions_b,
    #geometry data
    area, norm, dist,
    #other props 
    rho_a, rho_b, #kinda a u value because it changes with time but not explicitly tracked through u values
    species_molecular_weights
    )

    upwinded_mass_fractions = upwind(species_mass_fractions_a, species_mass_fractions_b, face_m_dot)

    du_species_mass_fractions_a -= face_m_dot * upwinded_mass_fractions
    du_species_mass_fractions_b += face_m_dot * upwinded_mass_fractions
end

#use get_cell_cp from helper functions to get cp for cell a and b

function enthalpy_advection!(
    #mutated vars
    du_temp_a, du_temp_b,
    #u values
    face_m_dot,
    temp_a, temp_b,
    #geometry data
    area, norm, dist,
    #other props 
    cp_a, cp_b,
    species_molecular_weights
    )
    cp_upwinded = upwind(cp_a, cp_b, face_m_dot)

    temp_upwinded = upwind(temp_a, temp_b, face_m_dot)

    energy_flux = face_m_dot * cp_upwinded * temp_upwinded 

    du_temp_a -= energy_flux
    du_temp_b += energy_flux
end


