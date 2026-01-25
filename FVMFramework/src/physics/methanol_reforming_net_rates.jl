
struct MSRReaction <: AbstractReaction
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

    adsorption_A_vec::Vector{Float64}
    adsorption_dH_vec::Vector{Float64}
end

struct MDReaction <: AbstractReaction
    heat_of_reaction::Float64
    delta_gibbs_free_energy::Float64
    K_gibbs_free_ref_temp::Float64
    kf_A::Float64 # Pre-exponential factor for forward reaction
    kf_Ea::Float64 # Activation energy for forward reaction
    reactants::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    reactant_stoich_coeffs::Vector{Int} #(stoich_coeff_for_reactant_1, stoich_coeff_for_reactant_2, etc...)
    products::Vector{Int} #(chemical_id_1, chemical_id_2, etc...)
    product_stoich_coeffs::Vector{Int} #[stoich_coeff_for_product_1, stoich_coeff_for_product_2, etc...]
    all_stoich_coeffs::Vector{Int} #(stoich_for_chemical_id_1, stoiach_for_chemical_id_2, etc...) #-1 = reactant, 1 = product

    adsorption_A_vec::Vector{Float64}
    adsorption_dH_vec::Vector{Float64}
end

struct WGSReaction <: AbstractReaction
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

    adsorption_A_vec::Vector{Float64}
    adsorption_dH_vec::Vector{Float64}
end


function PAM_reaction(rxn, molar_concentrations, T)
    pp = molar_concentrations .* ((R_gas * T) * 1e-5)
    #we divide by 1e-5 because PAM parameters are typically in bar
    # methanol, water, carbon_monoxide, hydrogen, carbon_dioxide
    P_CH3OH = pp[1]
    P_H2O = pp[2]
    P_CO = pp[3]
    P_H2 = max(pp[4], 1e-9)
    P_CO2 = pp[5]

    K_CH3O = van_t_hoff(rxn.adsorption_A_vec[1], rxn.adsorption_dH_vec[1], T)
    K_HCOO = van_t_hoff(rxn.adsorption_A_vec[2], rxn.adsorption_dH_vec[2], T)
    K_OH = van_t_hoff(rxn.adsorption_A_vec[3], rxn.adsorption_dH_vec[3], T)

    term_CH3O = K_CH3O * (P_CH3OH / sqrt(P_H2))
    term_HCOO = K_HCOO * (P_CO2 * sqrt(P_H2))
    term_OH = K_OH * (P_H2O / sqrt(P_H2))

    DEN = 1.0 + term_CH3O + term_HCOO + term_OH

    return P_CH3OH, P_H2O, P_CO, P_H2, P_CO2,   K_CH3O, K_HCOO, K_OH,   DEN
end


function net_reaction_rate(rxn::MSRReaction, molar_concentrations, T, kf_A, kf_Ea, kr_A, kr_Ea)
    P_CH3OH, P_H2O, P_CO, P_H2, P_CO2, 
    K_CH3O, K_HCOO, K_OH, 
    DEN = PAM_reaction(rxn, molar_concentrations, T)

    k_MSR = rxn.kf_A * exp(-rxn.kf_Ea / (R_gas * T))
    K_eq_MSR = K_gibbs_free(rxn.K_gibbs_free_ref_temp, T, rxn.delta_gibbs_free_energy, rxn.heat_of_reaction)
    driving_force = 1.0 - ( (P_CO2 * P_H2^3) / (K_eq_MSR * P_CH3OH * P_H2O) )
    
    rate = (k_MSR * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2) #max just to ensure we don't divide by zero

    return rate
end

function net_reaction_rate(rxn::MDReaction, molar_concentrations, T, kf_A, kf_Ea, kr_A, kr_Ea)
    P_CH3OH, P_H2O, P_CO, P_H2, P_CO2, 
    K_CH3O, K_HCOO, K_OH, 
    DEN = PAM_reaction(rxn, molar_concentrations, T)

    k_MD = rxn.kf_A * exp(-rxn.kf_Ea / (R_gas * T))
    K_eq_MD = K_gibbs_free(rxn.K_gibbs_free_ref_temp, T, rxn.delta_gibbs_free_energy, rxn.heat_of_reaction)
    driving_force = 1.0 - ( (P_CO * P_H2^2) / (K_eq_MD * P_CH3OH) )
    
    rate = (k_MD * K_CH3O * (P_CH3OH / sqrt(P_H2)) * driving_force) / (DEN^2)

    return rate
end

function net_reaction_rate(rxn::WGSReaction, molar_concentrations, T, kf_A, kf_Ea, kr_A, kr_Ea)
    P_CH3OH, P_H2O, P_CO, P_H2, P_CO2, 
    K_CH3O, K_HCOO, K_OH, 
    DEN = PAM_reaction(rxn, molar_concentrations, T)

    k_WGS = rxn.kf_A * exp(-rxn.kf_Ea / (R_gas * T))
    K_eq_WGS = K_gibbs_free(rxn.K_gibbs_free_ref_temp, T, rxn.delta_gibbs_free_energy, rxn.heat_of_reaction)
    driving_force = 1.0 - ( (P_CO2 * P_H2) / (K_eq_WGS * P_CO * P_H2O) )
    
    rate = (k_WGS * K_OH * ( (P_CO * P_H2O) / sqrt(P_H2) ) * driving_force) / (DEN^2)

    return rate
end