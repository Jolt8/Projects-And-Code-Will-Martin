using Unitful
using Symbolics, Nemo, Groebner
using Clapeyron
using NLsolve

motive_gas_mass_flow_inp = 0.2u"kg/s" #based on typical air compressor mass flow
motive_inlet_pressure_inp = 700u"kPa"
motive_inlet_temperature_inp = 293u"K"

suction_mass_flow_inp = 0.1u"kg/s" #based on an ER of 0.5 and motive mass flow of 0.2kg/s

cp_inp = 1.005u"kJ/(kg*K)"
cv_inp = 0.718u"kJ/(kg*K)"

function calc_choked_throat_area(motive_mass_flow, motive_inlet_pressure, motive_inlet_temperature, cp, cv)
    γ = cp / cv
    R_s = cp - cv
    return motive_mass_flow / motive_inlet_pressure * sqrt((R_s * motive_inlet_temperature / γ)) * (((γ + 1) / 2)^((γ + 1) / (2 * (γ - 1))))
end

choked_throat_area_unconverted = calc_choked_throat_area(motive_gas_mass_flow_inp, motive_inlet_pressure_inp, motive_inlet_temperature_inp, cp, cv)

choked_throat_area_found = uconvert(u"m^2", choked_throat_area_unconverted) #returns 0.00012100861121028794 m^2

sqrt(choked_throat_area_found / pi) #choked section circle radius, returns 0.006206306249421255 m
sqrt(choked_throat_area_found / pi) * 2 #choked section circle diameter, returns 0.01241261249884251  m



#this returns required final nozzle diameter (at the end of the nozzle) based on the choked nozzle diameter
function area_mach_relation(choked_throat_area, mach, cp, cv)
    γ = cp / cv
     return choked_throat_area * (1 / mach) * (((2 / (γ + 1)) * (1 + ((γ - 1) / 2)* mach^2)))^((γ + 1) / (2 * (γ - 1)))
end

#below is for ending area after choke
final_nozzle_area_found = area_mach_relation(choked_throat_area_found, 1.92, cp_inp, cv_inp) #returns 0.00019126112271002698  m^2

sqrt(final_nozzle_area_found / pi) #final section circle radius, returns 0.007802583303061411 m
sqrt(final_nozzle_area_found / pi) * 2 #final section circle diameter, returns 0.015605166606122822 m


#below is for starting area before choke
safety_factor = 1.1
final_nozzle_area_found = area_mach_relation(choked_throat_area_found, 0.3, cp_inp, cv_inp) * safety_factor #returns 0.0002708932076380669 m^2

sqrt(final_nozzle_area_found / pi) #final section circle radius, returns 0.009285902545861391 m
sqrt(final_nozzle_area_found / pi) * 2 #final section circle diameter, returns 0.018571805091722782 m



function nozzle_exit_static_pressure(inlet_pressure, mach, cp, cv)
    γ = cp / cv
    return inlet_pressure * (1 + ((γ - 1) / 2)* mach^2) ^ (-γ/(γ - 1))
end

nozzle_exit_static_pressure(motive_inlet_pressure_inp, 1.92, cp_inp, cv_inp)
