
function net_reaction_rate(chemical_reaction, molar_concentrations, T, kf_A, kf_Ea, kr_A, kr_Ea)
    R_gas_kj = 0.008314
    kf = (kf_A * exp(-kf_Ea / (R_gas_kj * T)))
    kr = (kr_A * exp(-kr_Ea / (R_gas_kj * T)))

    forward_term = 1.0
    for (i, species_id) in enumerate(chemical_reaction.reactants)
        concentration = molar_concentrations[species_id]
        stoich_coeff = chemical_reaction.reactant_stoich_coeffs[i]
        forward_term *= concentration^stoich_coeff
    end
    
    reverse_term = 1.0
    for (i, species_id) in enumerate(chemical_reaction.products)
        concentration = molar_concentrations[species_id]
        stoich_coeff = chemical_reaction.product_stoich_coeffs[i]
        reverse_term *= concentration^stoich_coeff
    end

    net_reaction_rate = ((kf * forward_term) - (kr * reverse_term))

    return net_reaction_rate
end

function K_gibbs_free(T_ref, T_actual, ΔG_rxn_ref, ΔH_rxn_ref)
    K_ref = exp(-ΔG_rxn_ref / (8.314e-3 * T_ref)) #R is in kJ

    ln_K_ratio = (-ΔH_rxn_ref / 8.314e-3) * (1/T_actual - 1/T_ref)
    
    K_T = K_ref * exp(ln_K_ratio)
    
    return K_T 
end

function react_cell!(
    #mutated vars
    du_mass_fractions, du_temp, 
    #caches (also mutated)
    molar_concentrations_cache, net_rates_cache, 
    #u data
    species_mass_fractions, cell_temp, 
    #geometry data
    vol,
    #props
    rho,
    #other data
    species_molecular_weights, cell_chemical_reactions_vec 
    )

    for (species_id, species_mass_fraction) in enumerate(species_mass_fractions)
        molar_concentrations_cache[species_id] = (rho * species_mass_fraction) / species_molecular_weights[species_id]
    end

    for (reaction_id, reaction) in enumerate(cell_chemical_reactions_vec)
        kf_A = reaction.kf_A
        kf_Ea = reaction.kf_Ea

        #find reverse pre exponential_factor
        K_ref = K_gibbs_free(reaction.K_gibbs_free_ref_temp, cell_temp, reaction.delta_gibbs_free_energy, reaction.heat_of_reaction)

        kr_A = (kf_A / K_ref) * exp(-reaction.heat_of_reaction / (8.314e-3 * cell_temp))

        #find reverse Ea
        kr_Ea = kf_Ea - reaction.heat_of_reaction

        net_rates_cache[reaction_id] = net_reaction_rate(reaction, molar_concentrations_cache, cell_temp, kf_A, kf_Ea, kr_A, kr_Ea)
    end

    for (reaction_id, reaction) in enumerate(cell_chemical_reactions_vec) 
        for (species_id, species_molar_concentration) in enumerate(molar_concentrations_cache)
            change_in_species_molar_concentration = 0.0

            for (reaction_idx, reaction) in enumerate(cell_chemical_reactions_vec)
                stoich = reaction.all_stoich_coeffs[species_id]
                change_in_species_molar_concentration += net_rates_cache[reaction_idx] * stoich

                du_temp[1] += net_rates_cache[reaction_idx] * -reaction.heat_of_reaction * vol
                #add heat of reaction to cell temp, - because -1.0 = exothermic
                #rate is in mol/m^3
                #du.temp[cell_id] is viewed as a 1 dimensional array
            end

            du_mass_fractions[species_id] += (change_in_species_molar_concentration * species_molecular_weights[species_id]) / rho
        end
    end
end
