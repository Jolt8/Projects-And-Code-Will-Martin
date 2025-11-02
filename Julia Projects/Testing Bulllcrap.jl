using Revise

using Unitful, DifferentialEquations
# using Clapeyron # Not used in this script
using Plots

using ModelingToolkit
using MethodOfLines
using DomainSets

# --- Constants (unchanged) ---
R = ustrip(uconvert(u"kJ/(mol*K)", 8.314u"J/(mol*K)"))
cross_sect_area = ustrip(uconvert(u"m^2", 0.3u"m^2"))
ρ_bulk = ustrip(uconvert(u"kg/m^3", 1000u"kg/m^3"))
Dp = ustrip(uconvert(u"m", 3u"mm"))
epsilon = 0.3 #Bed void fraction
mu_gas = ustrip(uconvert(u"bar*s", 2.0e-5u"Pa*s"))

MW_CH3OH = ustrip(uconvert(u"kg/mol", 0.03204u"g/mol"))
MW_H2O = ustrip(uconvert(u"kg/mol", 0.018015u"g/mol"))
MW_CO = ustrip(uconvert(u"kg/mol", 0.02801u"g/mol"))
MW_H2 = ustrip(uconvert(u"kg/mol", 0.002016u"g/mol"))
MW_CO2 = ustrip(uconvert(u"kg/mol", 0.04401u"g/mol"))

kf_A_SMR = ustrip(uconvert(u"mol/(bar^2*kg*s)", 300000u"mol/(bar^2*kg*s)"))
Ea_f_SMR = ustrip(uconvert(u"kJ/mol", 90u"kJ/mol"))
kr_A_SMR = ustrip(uconvert(u"mol/(bar^4*kg*s)", 0.008u"mol/(bar^4*kg*s)"))
Ea_r_SMR = ustrip(uconvert(u"kJ/mol", 43u"kJ/mol"))
Keq_ref_SMR = 1000
ΔH_SMR = ustrip(uconvert(u"kJ/mol", 49.4u"kJ/mol"))
ref_T_SMR = 298

kf_A_WGS = ustrip(uconvert(u"mol/(bar^2*kg*s)", 3000u"mol/(bar^2*kg*s)"))
Ea_f_WGS = ustrip(uconvert(u"kJ/mol", 60u"kJ/mol"))
kr_A_WGS = ustrip(uconvert(u"mol/(bar^2*kg*s)", 35000u"mol/(bar^2*kg*s)"))
Ea_r_WGS = ustrip(uconvert(u"kJ/mol", 100u"kJ/mol"))
Keq_ref_WGS= 36
ΔH_WGS = ustrip(uconvert(u"kJ/mol", -41.1u"kJ/mol"))
ref_T_WGS = 298

# --- PDE System Definition ---
@independent_variables t z 
@variables F_CH3OH(..) F_H2O(..) F_CO(..) F_H2(..) F_CO2(..) T(..) P(..)

Dt = Differential(t)
Dz = Differential(z)

# --- Intermediate Expressions (mostly unchanged) ---
F_total = F_CH3OH(t, z) + F_H2O(t, z) + F_CO(t, z) + F_H2(t, z) + F_CO2(t, z)

y_CH3OH = F_CH3OH(t, z) / F_total 
y_H2O = F_H2O(t, z) / F_total 
y_CO = F_CO(t, z) / F_total 
y_H2 = F_H2(t, z) / F_total 
y_CO2 = F_CO2(t, z) / F_total

p_CH3OH = y_CH3OH * P(t, z)
p_H2O = y_H2O * P(t, z)
p_CO = y_CO * P(t, z)
p_H2 = y_H2 * P(t, z)
p_CO2 = y_CO2 * P(t, z)

kf_SMR = kf_A_SMR * exp(-Ea_f_SMR / (R * T(t, z)))
Keq_SMR = Keq_ref_SMR * exp((-ΔH_SMR / R) * (1/T(t, z) - 1/ref_T_SMR))
kr_SMR = kf_SMR / Keq_SMR

kf_WGS = kf_A_WGS * exp(-Ea_f_WGS / (R * T(t, z)))
Keq_WGS = Keq_ref_WGS * exp((-ΔH_WGS / R) * (1/T(t, z) - 1/ref_T_WGS))
kr_WGS = kf_WGS / Keq_WGS

r_net_SMR = kf_SMR * (p_CH3OH * p_H2O) - kr_SMR * (p_CO * p_H2^3)
r_net_WGS = kf_WGS * (p_CO * p_H2O) - kr_WGS * (p_CO2 * p_H2)

# --- Define Reaction Source Terms ---
source_CH3OH = ρ_bulk * cross_sect_area * (-1 * r_net_SMR)
source_H2O = ρ_bulk * cross_sect_area * (-1 * r_net_SMR - 1 * r_net_WGS)
source_CO = ρ_bulk * cross_sect_area * (1 * r_net_SMR - 1 * r_net_WGS)
source_H2 = ρ_bulk * cross_sect_area * (3 * r_net_SMR + 1 * r_net_WGS)
source_CO2 = ρ_bulk * cross_sect_area * (1 * r_net_WGS)

# --- NEW: Define Interstitial Velocity ---
# This term converts the equation from a steady-state ODE in z to a transient PDE in t and z
interstitial_velocity = (R * T(t, z) * F_total) / (P(t, z) * epsilon * cross_sect_area)

# --- MODIFIED: Transient Mass Balance Equations ---
# Form: Dt(F) ~ velocity * (-Dz(F) + source)
eq_F_CH3OH = Dt(F_CH3OH(t, z)) ~ interstitial_velocity * (-Dz(F_CH3OH(t, z)) + source_CH3OH)
eq_F_H2O = Dt(F_H2O(t, z)) ~ interstitial_velocity * (-Dz(F_H2O(t, z)) + source_H2O)
eq_F_CO = Dt(F_CO(t, z)) ~ interstitial_velocity * (-Dz(F_CO(t, z)) + source_CO)
eq_F_H2 = Dt(F_H2(t, z)) ~ interstitial_velocity * (-Dz(F_H2(t, z)) + source_H2)
eq_F_CO2 = Dt(F_CO2(t, z)) ~ interstitial_velocity * (-Dz(F_CO2(t, z)) + source_CO2)

# --- Pressure Equation (unchanged, treated as steady-state) ---
MW_avg = y_CH3OH * MW_CH3OH + y_H2O * MW_H2O + y_CO * MW_CO + y_H2 * MW_H2 + y_CO2 * MW_CO2
mass_flow_total = F_CH3OH(t, z) * MW_CH3OH + F_H2O(t, z) * MW_H2O + F_CO(t, z) * MW_CO + F_H2(t, z) * MW_H2 + F_CO2(t, z) * MW_CO2
G = mass_flow_total / cross_sect_area
rho_gas = (P(t, z) * MW_avg) / (R * T(t, z))
ergun_term_1 = 150 * mu_gas * (1 - epsilon)^2 / (Dp^2 * epsilon^3)
ergun_term_2 = 1.75 * G * (1 - epsilon) / (Dp * epsilon^3)
ergun_expression = - (G / rho_gas) * (ergun_term_1 + ergun_term_2)
eq_P = Dz(P(t, z)) ~ ergun_expression

# --- MODIFIED: Temperature Equation ---
# Dt(T) ~ 0 represents an isothermal reactor where T at each point z is constant over time.
# A full energy balance would replace the 0 with heat of reaction and heat transfer terms.
eq_T = Dt(T(t, z)) ~ 0

eqs = [eq_F_CH3OH, eq_F_H2O, eq_F_CO, eq_F_H2, eq_F_CO2, eq_P, eq_T]

# --- Domains and Boundary/Initial Conditions (BCs/ICs) ---
t_start = 0.0
t_end = 200.0 # Reduced for faster initial testing

z_min = 0.0
z_max = 5.0

domains = [t ∈ Interval(t_start, t_end),
           z ∈ Interval(z_min, z_max)]

# BCs now clearly separate into initial conditions (at t=t_start)
# and boundary conditions (at z=z_min)
bcs = [
    # Initial conditions (state of the reactor at t=0)
    # Let's assume the reactor is initially filled with an inert gas or has very low reactant concentrations.
    F_CH3OH(t_start, z) ~ 1.0e-8,
    F_H2O(t_start, z) ~ 1.0e-8,
    F_CO(t_start, z) ~ 1.0e-8,
    F_H2(t_start, z) ~ 1.0e-8,
    F_CO2(t_start, z) ~ 1.0e-8,
    P(t_start, z) ~ 1.0,
    T(t_start, z) ~ 573.0,

    # Boundary conditions (feed to the reactor at z=0 for all t > 0)
    # Note: A step change in feed can cause stiffness. Let's use a smooth ramp-up.
    # ramp_factor = min(1.0, t/10.0 + 0.01) # Ramp up over 10 seconds
    F_CH3OH(t, z_min) ~ 0.1, # Using a constant feed for simplicity first
    F_H2O(t, z_min) ~ 0.1,
    F_CO(t, z_min) ~ 1.0e-8,
    F_H2(t, z_min) ~ 1.0e-8,
    F_CO2(t, z_min) ~ 1.0e-8,
    P(t, z_min) ~ 1.0,
    T(t, z_min) ~ 573.0,
]

@named pdesystem = PDESystem(eqs, bcs, domains, [t, z], [F_CH3OH(t, z), F_H2O(t, z), F_CO(t, z), F_H2(t, z), F_CO2(t, z), T(t, z), P(t, z)])

# --- Discretization and Solving ---
# Discretize only the spatial variable 'z'. Time is handled by the ODE solver.
# Increased number of points for better spatial resolution.
discretization = MOLFiniteDifference([z => 10], t) 
@time prob = discretize(pdesystem, discretization)

# Using a stiff solver is recommended for reactor models.
sol = solve(prob, Rodas5P(), saveat=1.0) # Using Rodas5P, a good stiff solver. saveat saves solution at specified time points.

# --- Plotting Results ---
z_grid = sol[z]
t_grid = sol[t]

# Extract solution for a specific variable, e.g., H2
F_H2_sol = sol[F_H2(t, z)]

# Plot 1: Concentration profile along the reactor at different times
plt_profiles = plot(z_grid, F_H2_sol[1, :], xlabel="Reactor Length (m)", ylabel="H2 Molar Flow (mol/s)", label="t = $(t_grid[1]) s", lw=2)
plot!(plt_profiles, z_grid, F_H2_sol[round(Int, end/4), :], label="t = $(t_grid[round(Int, end/4)]) s")
plot!(plt_profiles, z_grid, F_H2_sol[round(Int, end/2), :], label="t = $(t_grid[round(Int, end/2)]) s")
plot!(plt_profiles, z_grid, F_H2_sol[end, :], label="t = $(t_grid[end]) s (Final)")
title!("H2 Flow Profiles at Different Times")
display(plt_profiles)

# Plot 2: Outlet concentration over time (breakthrough curve)
outlet_z_index = last(eachindex(z_grid))
plt_breakthrough = plot(t_grid, F_H2_sol[:, outlet_z_index], xlabel="Time (s)", ylabel="Outlet H2 Molar Flow (mol/s)", label="H2 at Outlet", lw=2)
title!("Outlet H2 Flow vs. Time")
display(plt_breakthrough)

# Plot 3: 2D surface plot to see the full transient behavior
plt_surface = surface(z_grid, t_grid, F_H2_sol, xlabel="z (m)", ylabel="t (s)", zlabel="F_H2 (mol/s)", title="Surface Plot of H2 Molar Flow")
display(plt_surface)