using XLSX
using Unitful

xf = XLSX.readxlsx("C://Users//wille//Desktop//Legacy-Projects-And-Code-Will-Martin//Excel Projects//esterification_data_processing_for_julia.xlsx")

process_xf = xf[1]

timestamps = float.(vec(process_xf["A2:A10"])) .* u"s"

temperatures_deg_C = float.(vec(process_xf["B2:B10"])) .* u"°C"

temperature = uconvert.(u"K", temperatures_deg_C)

acetic_acid_moles = float.(vec(process_xf["D2:D10"])) .* u"mol"
ethanol_moles = float.(vec(process_xf["G2:G10"])) .* u"mol"
ethyl_acetate_moles = float.(vec(process_xf["I2:I10"])) .* u"mol"
water_moles = float.(vec(process_xf["H2:H10"])) .* u"mol"

function get_mass_fraction(species_moles, species_molecular_weight, rho, mixture_volume)
    (((species_moles / mixture_volume) * species_molecular_weight) / rho) 
end

mixture_volume = 50u"ml"
mixture_rho = 1000u"kg/m^3"

acetic_acid_mw = 60.05u"g/mol"
ethanol_mw = 46.07u"g/mol"
ethyl_acetate_mw = 88.11u"g/mol"
water_mw = 18.02u"g/mol"

acetic_acid_mass_fractions = get_mass_fraction.(acetic_acid_moles, acetic_acid_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
ethanol_mass_fractions = get_mass_fraction.(ethanol_moles, ethanol_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
ethyl_acetate_mass_fractions = get_mass_fraction.(ethyl_acetate_moles, ethyl_acetate_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
water_mass_fractions = get_mass_fraction.(water_moles, water_mw, mixture_rho, mixture_volume) .|> u"kg/kg"
