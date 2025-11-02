using Unitful
using Symbolics, Nemo, Groebner

@variables volumetric_flow, rotational_speed, rotor_width, eccentricity, rotor_diameter, stator_diameter

eqs = ((volumetric_flow / rotational_speed) ~ 2 * pi * eccentricity * (stator_diameter / 2) * rotor_width)

solved_for_variable = stator_diameter

sol_expression = symbolic_solve(eqs, solved_for_variable)[1]

#values to be subbed in
volumetric_flow_val = 0.000505u"m^3/s"
volumetric_flow_val = uconvert(u"m^3/s", volumetric_flow_val)

rotational_speed_rpm = 3000u"rpm"
rotational_speed_val= uconvert(u"s^-1", rotational_speed_rpm)

eccentricity_val = 1.5u"mm"
eccentricity_val = uconvert(u"m", eccentricity_val)

stator_diameter_val = 20u"mm"
#stator_diameter_val = uconvert(u"m^3/s", stator_diameter_val)

rotor_width_val = 20u"mm"
rotor_width_val = uconvert(u"m", rotor_width_val)


sub_dict = Dict(
    volumetric_flow => volumetric_flow_val,
    rotational_speed => rotational_speed_val,
    eccentricity => eccentricity_val,
    stator_diameter => stator_diameter_val,
    rotor_width => rotor_width_val
)

final_solution = substitute(sol_expression, sub_dict)
uconvert(u"mm", final_solution)


#NOTE: for future reference, if the final equation involves a √, you CAN'T use unitful units (mostly because it messes with area which is ^3 and unable to be square rooted)





