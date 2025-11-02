module Chemeqs
    export arrenhius_equation_pre_exponential_factor, arrenhius_equation_rate_constant, ChemicalReaction, HeterogeneousChemicalReaction, overall_rate, rate_law_for_chemicals, heterogeneous_rate_law_for_chemicals, K_eq, van_t_hoft, K_gibbs_free
    using Clapeyron
    using Unitful
    function arrenhius_equation_pre_exponential_factor(k, Ea, T)
        R_GAS = 8.314u"J/(mol*K)"
        result = (k / exp(-Ea / (R_GAS * T)))
        #result = ustrip(result)
        #return result * Unitful.unit(k)
        return result
    end

    function arrenhius_equation_rate_constant(A, Ea, T)
        R_GAS = 8.314u"J/(mol*K)"
        return (A * exp(-Ea / (R_GAS * T)))
    end

    struct ChemicalReaction
        name::String
        # Use Dicts for reactants/products with species name => stoichiometric coefficient
        # For rate law, these are the powers in elementary reactions.
        # For mass balance, they are the coefficients.
        reactants::Dict{String, Int}
        products::Dict{String, Int}
        #To maintain order throughout the process
        all_chemicals::Vector{String}
        kf_A::typeof(1.0u"L/(mol*s)") # Pre-exponential factor for forward reaction
        kf_Ea::typeof(1.0u"J/mol")   # Activation energy for forward reaction
        kr_A::typeof(1.0u"L/(mol*s)") # Pre-exponential factor for reverse reaction
        kr_Ea::typeof(1.0u"J/mol")   # Activation energy for reverse reaction
    end
    struct HeterogeneousChemicalReaction
        name::String
        # Use Dicts for reactants/products with species name => stoichiometric coefficient
        # For rate law, these are the powers in elementary reactions.
        # For mass balance, they are the coefficients.
        reactants::Dict{String, Int}
        products::Dict{String, Int}
        #To maintain order throughout the process
        all_chemicals::Vector{String}
        kf_A::typeof(1.0u"mol/(kg*s*Pa") # Pre-exponential factor for forward reaction
        kf_Ea::typeof(1.0u"J/mol")   # Activation energy for forward reaction
        kr_A::typeof(1.0u"mol/(kg*s*Pa)") # Pre-exponential factor for reverse reaction
        kr_Ea::typeof(1.0u"J/mol")   # Activation energy for reverse reaction
    end

    """Example for Steam Propane Reforming
    SPR_reaction = ChemicalReaction(
        "SPR",
        Dict("C3H8" => 1, "H2O" => 3),
        Dict("CO" => 3, "H2" => 7),
        ["C3H8", "H2O", "CO", "H2"],
        kf_A, Ea_f,
        kr_A, Ea_r
    )
    end


    """

    function overall_rate(reaction, model, temperature, pressure, concentrations_vec)
        total_moles = sum(concentrations_vec)
        mole_fractions = ustrip.(concentrations_vec ./ total_moles)

        activity_coeffs = activity_coefficient(model, ustrip(uconvert(u"kPa", pressure)), ustrip(uconvert(u"K", temperature)), mole_fractions)
        activities = activity_coeffs .* concentrations_vec
        
        reactants_activities = activities[1:length(reaction.reactants)]
        products_activities = activities[length(reaction.products)+1:end]

        forward_k = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
        reverse_k = arrenhius_equation_rate_constant(reaction.kr_A, reaction.kr_Ea, temperature)
        
        if typeof(reactants_activities) == Unitful.Quantity
            store_unit = unit(reactants_activities)
            reactants_activities = ustrip.(reactants_activities)
            products_activities = ustrip.(products_activities)
            return store_unit * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
        else
            return forward_k * prod(reactants_activities) - reverse_k * prod(products_activities)
        end
    end

    function rate_law_for_chemicals(reaction, model, temperature, pressure, concentrations_vec)
        total_moles = sum(concentrations_vec)
        mole_fractions = ustrip.(concentrations_vec ./ total_moles)

        activity_coeffs = activity_coefficient(model, ustrip(uconvert(u"kPa", pressure)), ustrip(uconvert(u"K", temperature)), mole_fractions)
        activities = activity_coeffs .* concentrations_vec

        reactants_activities = activities[1:length(reaction.reactants)]
        products_activities = activities[length(reaction.products)+1:end]

        forward_k = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
        reverse_k = arrenhius_equation_rate_constant(reaction.kr_A, reaction.kr_Ea, temperature)
        
        stoichs_and_signs = Dict()
        for chemical in reaction.all_chemicals
            if haskey(reaction.reactants, chemical)
                stoichs_and_signs[chemical] = -1 * reaction.reactants[chemical]
            else
                stoichs_and_signs[chemical] = 1 * reaction.products[chemical]
            end
        end

        rates = Dict()
        if typeof(reactants_activities) == Unitful.Quantity
            store_unit = unit(reactants_activities)
            reactants_activities = ustrip.(reactants_activities)
            products_activities = ustrip.(products_activities)

            for chemical in reaction.all_chemicals
                if haskey(reaction.reactants, chemical)
                    rates[chemical] = (store_unit * (stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))))
                else 
                    rates[chemical] = (store_unit * (stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))))
                end
            end
        else
            for chemical in reaction.all_chemicals
                if haskey(reaction.reactants, chemical)
                    rates[chemical] = stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
                else
                    rates[chemical] = stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
                end
            end
        end
        return rates
    end
    function heterogeneous_rate_law_for_chemicals(reaction, model, temperature, pressure, molar_flows_vec)
        total_molar_flow = sum(molar_flows_vec)
        mole_fractions = ustrip.(molar_flows_vec ./ total_molar_flow)

        activity_coeffs = activity_coefficient(model, ustrip(uconvert(u"kPa", pressure)), ustrip(uconvert(u"K", temperature)), mole_fractions)
        activities = activity_coeffs .* mole_fractions

        reactants_activities = activities[1:length(reaction.reactants)]
        products_activities = activities[length(reaction.products)+1:end]

        forward_k = arrenhius_equation_rate_constant(reaction.kf_A, reaction.kf_Ea, temperature)
        reverse_k = arrenhius_equation_rate_constant(reaction.kr_A, reaction.kr_Ea, temperature)
        
        stoichs_and_signs = Dict()
        for chemical in reaction.all_chemicals
            if haskey(reaction.reactants, chemical)
                stoichs_and_signs[chemical] = -1 * reaction.reactants[chemical]
            else
                stoichs_and_signs[chemical] = 1 * reaction.products[chemical]
            end
        end

        rates = Dict()
        if typeof(reactants_activities) == Unitful.Quantity
            store_unit = unit(reactants_activities)
            reactants_activities = ustrip.(reactants_activities)
            products_activities = ustrip.(products_activities)

            for chemical in reaction.all_chemicals
                if haskey(reaction.reactants, chemical)
                    rates[chemical] = (store_unit * (stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))))
                else 
                    rates[chemical] = (store_unit * (stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))))
                end
            end
        else
            for chemical in reaction.all_chemicals
                if haskey(reaction.reactants, chemical)
                    rates[chemical] = stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
                else
                    rates[chemical] = stoichs_and_signs[chemical] * (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
                end
            end
        end
        return rates
    end

    function K_eq(forward_k, reverse_k)
        return forward_k / reverse_k
    end

    #K1 is initial equilibrum constant, returns equilibrium constant after temperature change
    function van_t_hoft(K_eq1, T1, T2, heat_of_reaction)
        R_GAS = 8.314u"J/(mol*K)"
        return K_eq1 * exp((-heat_of_reaction / R_GAS) * (1/T2 - 1/T1))
    end
    
    function Q_gibbs_free(reaction, clapeyron_model, T, p, z, heat_of_reaction)
        potentials = chemical_potential(clapeyron_model, p, T, z)
        R_gas = 8.31446261815324u"J/(K*mol)"

        delta_gibbs = 0u"kJ"
        i = 0
        for chemical in reaction.all_chemicals
            if haskey(reaction.reactants, chemical) 
                delta_gibbs += potentials[i]
            else 
                delta_gibbs -= potentials[i]
            i += 1
            end
        end
        return exp((-1*(reactant_gibbs_free_energy - product_gibbs_free_energy)) / (R_gas * T))
        #return exp((heat_of_reaction - (T * (reactant_entropy - product_entropy))) / (R_gas * T))
    end  
    function K_gibbs_free(T_ref, T_actual, ΔG_rxn_ref, ΔH_rxn_ref)
        T_actual = uconvert(u"K", T_actual) # Ensure T is a plain number in Kelvin

        R = 8.314e-3u"kJ/(mol*K)" # Gas constant in kJ
        
        K_ref = exp(-ΔG_rxn_ref / (R * T_ref))

        ln_K_ratio = (-ΔH_rxn_ref / R) * (1/T_actual - 1/T_ref)
        
        K_T = K_ref * exp(ln_K_ratio)
        
        return ustrip(K_T) # Return a unitless equilibrium constant
    end
end # module Chemeqs

