const R_gas = 8.314

function upwind(var_a, var_b, mass_flow_rate)
    if mass_flow_rate > 0.0
        return var_a
    else 
        return var_b
    end
end

function harmonic_mean(val_a, val_b)
    2 * val_a * val_b / (val_a + val_b)
end

function get_cell_rho(
    mw_avg_cache_for_cell, #cached vars
    pressure, temp, species_mass_fractions, #u values
    species_molecular_weights #other props 
    )
    mw_avg_cache_for_cell = 0.0

    for i in eachindex(species_molecular_weights)
        mw_avg_cache_for_cell += species_mass_fractions_a[i] / species_molecular_weights[i]
    end

    mw_avg_cache_for_cell = mw_avg_cache_for_cell^-1.0

    return (pressure * avg_mw_a) / (R_gas * temp_a)
end

function get_cell_cp(
        cp_avg_cache_for_cell, #cached vars
        species_mass_fractions, #u values
        species_cps #other props
    )
    cp_avg_cache_for_cell = 0.0

    for i in eachindex(species_mass_fractions)
        cp_avg_cache_for_cell += species_mass_fractions[i] * species_cps[i]
    end

    return cp_avg_cache_for_cell
end
