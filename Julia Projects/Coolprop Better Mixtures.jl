using CoolProp

species = ["Water", "Ethanol"]

mole_fractions = [0.5, 0.5]

function coolprop_mixture(species, mole_fractions)
    total_fraction = sum(mole_fractions)
    corrected_mole_fractions = [mole_fraction / total_fraction for mole_fraction in mole_fractions]

    coolprop_mixture_string = ""
    for i in eachindex(species)
        if i != length(species)
            coolprop_mixture_string *= species[i] * "[" * string(corrected_mole_fractions[i]) * "]" * "&"
        else 
            coolprop_mixture_string *= species[i] * "[" * string(corrected_mole_fractions[i]) * "]"
        end
    end
    return coolprop_mixture_string
end

coolprop_mixture(species, mole_fractions)