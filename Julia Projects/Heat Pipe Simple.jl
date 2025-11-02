using Unitful, DifferentialEquations
using Plots

using ModelingToolkit
using MethodOfLines
using DomainSets

#using Clapeyron

Cp = ustrip(uconvert(u"J/(g*K)", 0.385u"J/(g*K)"))
P = ustrip(uconvert(u"Pa", 1u"bar"))

k = ustrip(uconvert(u"W/(m*K)", 401u"W/(m*K)"))

@independent_variables t x
@variables T_j(..) 
    
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

rho = ustrip(uconvert(u"g/m^3", 8960u"kg/m^3"))
#P / (461.5 * T_j(t, x))
heat_input = 10000.0
T_eq = rho * Cp * Dt(T_j(t, x)) ~ k * Dxx(T_j(t, x)) + heat_input

t_start = 0
t_end = 2000000.0

x_min = 0
x_max = 1

eqs = [T_eq]

domains = [t ∈ Interval(t_start, t_end),
           x ∈ Interval(x_min, x_max)]
bcs = [
    T_j(t_start, x) ~ 300,
    T_j(t, x_min) ~ 300,
]

@named pdesystem = PDESystem(eqs, bcs, domains, [t, x], [T_j(t, x)])

x_points = 6 #This is really weird behaviour but setting this to 5 makes the simulaton fail to run because it says BoundsError: attempt to access 5-element Vector{SymbolicUtils.BasicSymbolic{Real}} at index [CartesianIndex{1}[CartesianIndex(0,), CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,), CartesianIndex(4,), CartesianIndex(5,)]]
discretization = MOLFiniteDifference([x => x_points], t) # Discretize z into 20 points
@time prob = discretize(pdesystem, discretization)

sol = solve(prob, Rosenbrock23())
x_grid = sol[x]
t_grid = sol[t]
T_j_sol = sol[T_j(t, x)]

desired_max_time = 20000
max_time_idx = argmin(abs.(t_grid .- desired_max_time))
time_at_max = t_grid[max_time_idx] #* 1.0u"s"
plot(t_grid[1:max_time_idx], sol[T_j(t, x)][1:max_time_idx, 2], xlabel="time (s)", label=false, ylabel="Temperature (K)", lw=2)
sol[T_j(t, x)][1:max_time_idx, 2]

time_for_plot = time_at_max
current_time = t_grid[max_time_idx]
x_s = 1
x_e = length(x_grid)
x_slice = x_grid[x_s:x_e]
plot(x_slice, sol[T_j(t, x)][max_time_idx, x_s:x_e], xlabel="Heat Pipe Length (m)", label=false, ylabel="Temperature (K)", lw=2)
sol[T_j(t, x)][max_time_idx, x_s:x_e]
