
struct PowerLawReaction <: AbstractReaction
    heat_of_reaction::Float64
    delta_gibbs_free_energy::Float64
    K_gibbs_free_ref_temp::Float64
    kf_A::Float64 # Pre-exponential factor for forward reaction
    kf_Ea::Float64 # Activation energy for forward reaction
    reactants::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    reactant_stoich_coeffs::Vector{Int} #(stoich_coeff_for_reactant_1, stoich_coeff_for_reactant_2, etc...)
    products::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    product_stoich_coeffs::Vector{Int} #[stoich_coeff_for_product_1, stoich_coeff_for_product_2, etc...]
    all_stoich_coeffs::Vector{Int} #(stoich_for_chemical_id_1, stoich_for_chemical_id_2, etc...) #-1 = reactant, 1 = product
end

function net_reaction_rate(chemical_reaction::PowerLawReaction, molar_concentrations, T, kf_A, kf_Ea, kr_A, kr_Ea)
    kf = (kf_A * exp(-kf_Ea / (R_gas * T)))
    kr = (kr_A * exp(-kr_Ea / (R_gas * T)))

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
    K_ref = exp(-ΔG_rxn_ref / (R_gas * T_ref))

    ln_K_ratio = (-ΔH_rxn_ref / R_gas) * (1/T_actual - 1/T_ref)
    
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
    species_molecular_weights, reactions, cell_kg_cat_per_m3_for_each_reaction
    )

    for (species_id, species_mass_fraction) in enumerate(species_mass_fractions)
        molar_concentrations_cache[species_id] = (rho * species_mass_fraction) / species_molecular_weights[species_id]
    end

    for (reaction_id, reaction) in enumerate(reactions)
        kf_A = reaction.kf_A
        kf_Ea = reaction.kf_Ea

        #find reverse pre exponential_factor
        K_ref = K_gibbs_free(reaction.K_gibbs_free_ref_temp, cell_temp, reaction.delta_gibbs_free_energy, reaction.heat_of_reaction)

        kr_A = (kf_A / K_ref) * exp(-reaction.heat_of_reaction / (R_gas * cell_temp)) 

        #find reverse Ea
        kr_Ea = kf_Ea - reaction.heat_of_reaction

        net_rates_cache[reaction_id] = net_reaction_rate(reaction, molar_concentrations_cache, cell_temp, kf_A, kf_Ea, kr_A, kr_Ea) * cell_kg_cat_per_m3_for_each_reaction[reaction_id]
        #rate returned by net_reactions_rate is in mol / (kg_cat * s), so we times by the cell's kg_cat / m3
        #rate is now in mol / (m3 * s)
    end
    
    for (species_id, species_molar_concentration) in enumerate(molar_concentrations_cache)
        change_in_species_molar_concentration = 0.0

        for (reaction_id, reaction) in enumerate(reactions)
            stoich = reaction.all_stoich_coeffs[species_id]
            change_in_species_molar_concentration += net_rates_cache[reaction_id] * stoich
        end

        du_mass_fractions[species_id] += (change_in_species_molar_concentration * species_molecular_weights[species_id]) / rho
        # rate (mol/(m3*s)) * MW (g/mol) / rho (g/m3) = unitless/s
    end

    for (reaction_id, reaction) in enumerate(reactions) 
        du_temp[1] += net_rates_cache[reaction_id] * (-reaction.heat_of_reaction) * vol
        # rate (mol/(m3*s)) * H (J/mol) * vol (m3) = J/s = Watts
    end
end