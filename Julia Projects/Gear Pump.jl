using Unitful
using Symbolics, Nemo, Groebner


@variables number_of_teeth, gear_width_to_outer_diameter_ratio, volumetric_flow, rotational_speed, gear_modulus, gear_width, pitch_diameter, outer_diameter

eqs = [
    outer_diameter ~ gear_modulus * (2 + number_of_teeth),
    pitch_diameter ~ gear_modulus * number_of_teeth,
    (volumetric_flow / rotational_speed) ~ (pi / 2) * gear_width * (outer_diameter^2 - pitch_diameter^2)
]

sub_map = Dict(
    outer_diameter => eqs[1].rhs,
    pitch_diameter => eqs[2].rhs
)

expanded_volumetric_flow_eq = substitute(eqs[3], sub_map)


solved_for_variable = gear_width


sol_expression = symbolic_solve(expanded_volumetric_flow_eq, solved_for_variable)[1]

#gear_width_to_outer_diameter_ratio_val = 1.5

volumetric_flow_val = 0.0023u"m^3/s"

gear_width_val = 1.5u"mm"

rotational_speed_rpm = 300u"rpm"
rotational_speed_val= uconvert(u"s^-1", rotational_speed_rpm)
#rotational_speed_val = 300

gear_modulus_mm = 1.5u"mm"
gear_modulus_val = uconvert(u"m", gear_modulus_mm)

number_of_teeth_val = 12

pitch_diameter_val = gear_modulus_val * number_of_teeth_val
outer_diameter_val = gear_modulus_val * (2 + number_of_teeth_val)


sub_dict = Dict(
    volumetric_flow => volumetric_flow_val,
    rotational_speed => rotational_speed_val,
    gear_modulus => gear_modulus_val,
    number_of_teeth => number_of_teeth_val,
    gear_width => gear_width_val,
)


final_solution = substitute(sol_expression, sub_dict)

#NOTE: for future reference, if the final equation involves a √, you CAN'T use unitful units





