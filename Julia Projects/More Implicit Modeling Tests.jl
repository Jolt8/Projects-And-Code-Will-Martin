using ModelingToolkit
using Symbolics, Nemo, Groebner
using Unitful
using MacroTools


#BELOW IS ONLY TOOLING
#thin walled 
@variables wall_thickness, pressure, pipe_inner_diameter, material_yield_strength, joint_efficiency_factor
eq = wall_thickness ~ (pressure * pipe_inner_diameter) / (2 * joint_efficiency_factor * (8//10 * material_yield_strength))

sol_expression = symbolic_solve(eq, wall_thickness)[1]
#this gives (5pipe_inner_diameter*pressure) / (8joint_efficiency_factor*material_yield_strength which is not very readable

eq = wall_thickness ~ (pressure * pipe_inner_diameter) / (2 * joint_efficiency_factor * (0.8 * material_yield_strength))
sol_expression = solve_for(eq, wall_thickness, simplify=false)
#this gives (pipe_inner_diameter*pressure) / (1.6joint_efficiency_factor*material_yield_strength)
#for some reason solving for material_yield_strength errors with AssertionError: islinear


#Below will be the only things evaluated at runtime
function get_wall_thickness(pressure, pipe_inner_diameter, joint_efficiency_factor, material_yield_strength)
    return (pipe_inner_diameter*pressure) / (1.6joint_efficiency_factor*material_yield_strength)
end

#Heat Duty Energy Balance
@variables heat_load, mass_flow_h, specific_heat_capacity_h, T_in_h, T_out_h, T_in_c, T_out_c, mass_flow_c, specific_heat_capacity_c, T_out_c, T_in_c
begin eq = 
    heat_load ~ 
    mass_flow_h * specific_heat_capacity_h * (T_in_h - T_out_h) ~
    mass_flow_c * specific_heat_capacity_c * (T_out_c - T_in_c)
end

solve_for(eq, heat_load, simplify=false)

@variables U, A, F, test
#Requred area 
function get_LMTD(T_in_h, T_out_h, T_in_c, T_out_c)
    dT1 = T_in_h - T_out_c
    dT2 = T_out_h - T_in_c
    return (delta_T1 - delta_T2) / log(delta_T1 / delta_T2)
end

#dT1 = T_in_h - T_out_c
#dT2 = T_out_h - T_in_c
#test2 = (delta_T1 - delta_T2) / log(delta_T1 / delta_T2)
#not sure if the function or straight code is better, but they both work

begin eq = 
    heat_load ~ U * A * F * get_LMTD(T_in_h, T_out_h, T_in_c, T_out_C)
end

try 
    symbolic_linear_solve(eq, A)[1]
catch
    symbolic_solve(eq, A)[1]
end

