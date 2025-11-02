using Unitful
using Symbolics, Nemo, Groebner
using Clapeyron
using LinearAlgebra
using LinearSolve
using Roots

function fenske_underwood_gilliland_method(clapeyron_model, T, P, respective_mole_fractions_in_feed, feed_molar_flow, LK_mole_fraction_in_distillate, HK_mole_fraction_in_bottoms, actual_reflux_ratio)
    #note that if actual reflux ratio is none, it will return the minimum reflux ratio so you can build on that 
    Nc = length(respective_mole_fractions_in_feed)


    #get the saturation pressures of each component
    saturation_pressures = []
    for i in range(1, length(respective_mole_fractions_in_feed))
        current_model = PR([clapeyron_model.components[i]])
        push!(saturation_pressures, saturation_pressure(current_model, T)[1])
    end
    
    #get the names of each component
    component_names = clapeyron_model.components
    
    #zip them up in a dict so they can be sorted by saturation pressures and maintain their respective component
    zipped_vals = []
    for i in range(1, length(saturation_pressures))
        push!(zipped_vals, [saturation_pressures[i], component_names[i]])
    end

    #sort them
    sorted_pairs = sort(zipped_vals, by = x -> x[1])

    #unpack them from the dict into a new array
    for i in range(1, length(sorted_pairs))
        saturation_pressures[i] = sorted_pairs[i][1]
        component_names[i] = sorted_pairs[i][2]
    end

    light_key_p_sat = saturation_pressures[end] # Higher volatility
    heavy_key_p_sat = saturation_pressures[1]  # Lower volatility
    
    volatility_ratios = saturation_pressures ./ heavy_key_p_sat
    
    relative_volatility = light_key_p_sat / heavy_key_p_sat #usually denoted by the symmbol α
        
    HK_mole_fraction_in_distillate = 1 - HK_mole_fraction_in_bottoms
    LK_mole_fraction_in_bottoms = 1 - LK_mole_fraction_in_distillate
    
    #Which of these two is cleaner?
    N_min = (
        log10((LK_mole_fraction_in_distillate / HK_mole_fraction_in_distillate) * (HK_mole_fraction_in_bottoms / LK_mole_fraction_in_bottoms))
        / (log10(relative_volatility))
    )

    N_min = ((log10((LK_mole_fraction_in_distillate / (1 - LK_mole_fraction_in_distillate)) * (HK_mole_fraction_in_bottoms / (1 - HK_mole_fraction_in_bottoms))))
                / (log10(relative_volatility)))
    
    #minimum_reflux_ratio 

    
    A = log10((1 - HK_mole_fraction_in_bottoms) / HK_mole_fraction_in_bottoms)
    #calculate the recovery ratios for all other components
    distillate_recovery_ratios = []
    for ratio in volatility_ratios
        push!(distillate_recovery_ratios, (10^A) * (ratio^N_min) / (1 + (10^A) * (ratio^N_min)))
    end
    
    bottoms_recovery_ratios = 1 .- distillate_recovery_ratios
    
    respective_mass_flows = respective_mole_fractions_in_feed .* feed_molar_flow
        
    #to find underwood constant
    function find_underwood_constant(underwood_constant)
        results1 = []
        for i in range(1, length(Nc))
            push!(results1, ((volatility_ratios[i] * respective_mole_fractions_in_feed[i]) / (volatility_ratios[i] - underwood_constant)))
        end
        return sum(results1) + 1 - 0.5
    end
    max_retries1 = 500
    retries1 = 0
    found_underwood_constant = nothing
    while isnothing(found_underwood_constant) && retries1 < max_retries1
        found_underwood_constant = find_zero(find_underwood_constant, 1.5)
        retries1 += 1
    end
    
    #to find minimum reflux using underwood 
    function underwood_minimum_reflux(minimum_reflux_ratio)
        results2 = []
        for i in range(1, length(Nc))
            push!(results2, (volatility_ratios[i] * distillate_recovery_ratios[i]) / (volatility_ratios[i] - minimum_reflux_ratio))
        end
        return sum(results2) - minimum_reflux_ratio - 1
    end
    max_retries2 = 500
    retries2 = 0
    found_minimum_reflux_ratio = nothing
    while isnothing(found_minimum_reflux_ratio) && retries2 < max_retries2
        found_minimum_reflux_ratio = find_zero(underwood_minimum_reflux, 1.5)
        retries2 += 1
    end

    X = (actual_reflux_ratio - found_minimum_reflux_ratio) / (actual_reflux_ratio + 1)
    
    #Y = 1 - math.exp(((1 + (54.4*X))/(11 + (117.2*X)))*((X-1)/(X)))
    
    Y = (1 - (X^0.0031)) / (1 - (0.99357 * (X^0.0031)))
    
    #Gilland_correlation
    N_actual = (N_min + Y) / (1 - Y)
    
    
    if N_actual == 1
        println("")
        println("")
        println("Your inputted actual reflux ratio:", actual_reflux_ratio, "is less than the calculated:", found_minimum_reflux_ratio, "please increase it!")
        println("")
        println("")
    end
    
    println("compound_names: ", component_names, " saturation_pressures: ", saturation_pressures,)
    println("volatility_ratios: ", volatility_ratios,)
    println("distillate_recovery_ratios: ", distillate_recovery_ratios, " bottoms_recovery_ratios: ", bottoms_recovery_ratios,)      
    println("N_min: ", N_min, " underwood_constant: ", found_underwood_constant, " minimum reflux: ", found_minimum_reflux_ratio,)          
    println("light key psat: ", light_key_p_sat, " heavy_key psat: ", heavy_key_p_sat, )          
    println("min_reflux_ratio: ", found_minimum_reflux_ratio, " actual_reflux_ratio: ", actual_reflux_ratio, )         
    println("N_min: ", N_min, " N_actual: ", N_actual)
    println("")
    
    return("N_actual: ", N_actual, " min_reflux_ratio: ", found_minimum_reflux_ratio, " actual_reflux_ratio: ", actual_reflux_ratio)
end


model = PR(["Acetic Acid", "Ethanol", "Water", "Ethyl Acetate"])


fenske_underwood_gilliland_method(
    model, #clapeyron_model
    300, #T
    100000, #P
    [0.1, 1, 5, 5], #respective_mole_fractions_in_feed
    100, #feed_molar_flow
    0.60, #LK_mole_fraction_in_distillate=
    0.60, #HK_mole_fraction_in_bottoms=
    1 #actual_reflux_ratio=
)
