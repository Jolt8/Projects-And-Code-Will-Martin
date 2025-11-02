function rate_law_for_selected_chemical(selected_chemical, chemical_reaction, forward_k, reverse_k, activity_coeffs, concentrations, coefficients)
    dividing_signs = [
        "<>", "<->", "<-->", "⇋", "⇆", 
        "-->", "->", ">"
    ]
    reactants_side = nothing
    products_side = nothing
    for i in eachindex(dividing_signs)
        if occursin(dividing_signs[i], chemical_reaction) == true
            reactants_side = split(chemical_reaction, dividing_signs[i])[1]
            products_side = split(chemical_reaction, dividing_signs[i])[2]
            break
        end
    end
    rate_law_sign = 1

    num_of_reactants = count(==('+'), reactants_side) + 1
    num_of_products = count(==('+'), products_side) + 1

    selected_chemical_position = 0
    #if the chemical is a reactant, make it the opposite
    if occursin(selected_chemical, reactants_side)
        rate_law_sign *= -1
        
        target_range = findfirst("+", reactants_side)
        if isnothing(target_range)
            return 0
        end
        left_part = reactants_side[1 : first(target_range) - 1]
        chemical_position = count("+", left_part)
        selected_chemical_position = chemical_position + 1

    #if it's a product, keep it unchanged 
    else
        rate_law_sign *= 1
        target_range = findfirst("+", reactants_side)
        if isnothing(target_range)
            return 0
        end
        left_part = products_side[1 : first(target_range) - 1]
        chemical_position = count("+", left_part)
        selected_chemical_position = chemical_position + num_of_reactants + 1
        #not strictly necessary
    end

    activities = activity_coeffs .* concentrations

    reactants_activities = activities[1:num_of_reactants]
    products_activities = activities[num_of_reactants+1:end]
    
    rate_constant = rate_law_sign * coefficients[selected_chemical_position] *  (forward_k * prod(reactants_activities) - reverse_k * prod(products_activities))
    return rate_constant
end


selected_chemical = "Acetic Acid"

chemical_reaction = "Acetic Acid + Ethanol <--> Ethyl Acetate + Water"

activity_coefficients = activity_coefficient(model, 100000.0, 300.0, z)

concentrations = z .* u"mol/L" 

rate_law_for_selected_chemical(selected_chemical, chemical_reaction, 8.0e-5u"L/(mol*s)", 5.0e-5u"L/(mol*s)", activity_coefficients, concentrations, [1,1,1,1])