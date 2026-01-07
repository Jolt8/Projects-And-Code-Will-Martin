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