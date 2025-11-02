using Clapeyron
using Unitful
    
SMR_reaction = ChemicalReaction(
    "SMR",
    Dict("CH3OH" => 1, "H2O" => 1),
    Dict("CO" => 1, "H2" => 3),
    ["CH3OH", "H2O", "CO", "H2"],
    1.0u"L/(mol*s)", 1.0u"kJ/mol",
    1.0u"L/(mol*s)", 1.0u"kJ/mol"
)

clapeyron_model = PR(["Methanol", "water", "carbon monoxide", "hydrogen"])

reactants_model = PR(clapeyron_model.components[1:length(SMR_reaction.reactants)])

products_model = PR(clapeyron_model.components[length(SMR_reaction.reactants)+1:end])


T = 300u"°C" |> u"K"

p = 1u"atm" |> u"Pa"

z_all = [1, 1, 1, 1]


R_gas = 8.31446261815324u"J/(K*mol)"

potentials = chemical_potential(clapeyron_model, p, T, z)
R_gas = 8.31446261815324u"J/(K*mol)"

delta_gibbs = 0u"J/mol"

delta_gibbs += SMR_reaction.products["H2"] * potentials[4]
i = 1
for chemical in SMR_reaction.all_chemicals
    if haskey(SMR_reaction.reactants, chemical) 
        delta_gibbs += SMR_reaction.reactants[chemical] * potentials[i]
    else 
        delta_gibbs -= SMR_reaction.products[chemical] * potentials[i]
    i += 1
    end
end
    
exp((-1*(delta_gibbs)) / (R_gas * T))

(1/(300u"K") - 1/(300u"K"))



heat_of_reaction = 100u"kJ/mol"

#return exp((heat_of_reaction - (T * (reactant_entropy - product_entropy))) / (R_gas * T))
