using Revise

using Unitful, DifferentialEquations
using Plots

using ModelingToolkit
using MethodOfLines
using DomainSets

using Clapeyron
fluid_model = PR(["Sodium"])

evap_length = ustrip(uconvert(u"m", 72u"cm"))
adia_length = ustrip(uconvert(u"m", 10u"cm"))
cond_length = ustrip(uconvert(u"m", 100u"cm"))
pipe_length = sum([evap_length, adia_length, cond_length])
pipe_outer_radius = ustrip(uconvert(u"m", 3u"cm"))
wick_outer_radius = ustrip(uconvert(u"m", 2.8u"cm"))
wick_inner_radius = ustrip(uconvert(u"m", 2.7u"cm"))
pipe_tilt_angle = ustrip(uconvert(u"°", 0u"°"))

wick_cs_a = pi * (wick_outer_radius^2 - wick_inner_radius^2)
vapor_core_cs_a = pi * wick_inner_radius^2
wall_cs_a = pi * (pipe_outer_radius^2 - wick_outer_radius^2)

# Fluid represents sodium in our case
rho_l_func(T) = 1011.6 - 0.2208*T - 1.46e-5*T^2 + 5.69e-9*T^3 # Liquid density (kg/m^3)
rho_v_func(T, P) = P / (361.3 * T) # Vapor density (kg/m^3) using ideal gas law, R_gas = 361.3 J/kg-K

mu_l_func(T) = 2.415e-5 * 10^(248.6 / (T - 140.0)) # Liquid viscosity (Pa*s)
mu_v_func(T) = 1.002e-6 + 1.986e-8*T + 4.58e-12*T^2 # Vapor viscosity (Pa*s)

k_l_func(T) = 124.4 - 0.1138*T + 5.523e-5*T^2 - 1.184e-8*T^3 # Liquid conductivity (W/m-K)
H_vap_l_func(T) = (4.9886e6 - 5.751e2*T - 2.158*T^2) / 1000 # Latent heat (J/kg), converted from J/kmol

dynam_visc_l_func(T) = 0.2323 - 1.01e-4 * T # Surface tension (N/m)

# Wick and Wall properties
wall_k = ustrip(uconvert(u"W/(m*K)", 20.0u"W/(m*K)"))

saturated_wick_effective_k = ustrip(uconvert(u"W/(m*K)", 40u"W/(m*K)"))
wick_permeability = ustrip(uconvert(u"m^2", 1e-11u"m^2"))
wick_porosity = 0.6
g = ustrip(uconvert(u"m/s^2", 9.80665u"m/s^2"))
wick_effective_capillary_radius = ustrip(uconvert(u"m", 50u"μm")) #length usually in mm or μm

P_ref = ustrip(uconvert(u"Pa", 101325u"Pa"))
T_ref = ustrip(uconvert(u"K", 1156.0u"K"))
fluid_specific_gas_constant = ustrip(uconvert(u"J/(kg*K)", 361.3u"J/(kg*K)"))


@independent_variables z
@syms T_v(z)::Real P_v(z)::Real P_l(z)::Real m_dot_v(z)::Real m_dot_l(z)::Real q_axial(z)::Real m_dot_per_length(z)::Real #rho_l(..) rho_v(..) mu_l(..) mu_v(..) k_l(..) H_vap(..) dynam_visc_l(..)
#@parameters rho_l(..) rho_v(..) mu_l(..) mu_v(..) H_vap_l(..) Re_v(..)
Dz = Differential(z)

rho_l = rho_l_func(T_v(z))
rho_v = rho_v_func(T_v(z), P_v(z))
mu_l = mu_l_func(T_v(z))
mu_v = mu_v_func(T_v(z))
H_vap_l = H_vap_l_func(T_v(z))
dynam_visc_l = dynam_visc_l_func(T_v(z))

Re_v = abs(m_dot_v(z)) * (2 * wick_inner_radius) / (vapor_core_cs_a * mu_v)
# Fanning friction factor (using a smooth transition from laminar to turbulent)
f_v = 16
#ifelse(Re_v < 2300, 16 / Re_v, 0.079 / (Re_v^0.25))

thermal_res_wall = log(pipe_outer_radius / wick_outer_radius) / (2 * pi * wall_k)
thermal_res_wick = log(wick_outer_radius / wick_inner_radius) / (2 * pi * saturated_wick_effective_k)

Q_total_in = 1000.0 # Total power input in Watts
q_evap = Q_total_in / evap_length
q_cond = -Q_total_in / cond_length
q_input(z) = 1
#ifelse(z <= evap_length, q_evap,
                #ifelse(z <= evap_length + adia_length, 0.0, q_cond))

#Defines m_dot_per_length based on axial heat balance
#q_input(z) is external heat added/removed per unit length
eq_energy_balance = Dz(q_axial(z)) ~ q_input(z) - (m_dot_per_length(z) * H_vap_l)

# Axial heat conduction
eq_axial_conduction = q_axial(z) ~ - (wall_k * wall_cs_a + saturated_wick_effective_k * wick_cs_a) * Dz(T_v(z))

#Vapor mass conservation
eq_m_dot_v = Dz(m_dot_v(z)) ~ m_dot_per_length(z)

#Liquid mass conservation
eq_m_dot_l = Dz(m_dot_l(z)) ~ -m_dot_per_length(z)

friction_term = (2 * f_v / (2 * wick_inner_radius)) * m_dot_v(z) * abs(m_dot_v(z)) / (vapor_core_cs_a^2 * rho_v)

axial_momentum_flux_correc_factor = 1
#momentum_term = (axial_momentum_flux_correc_factor / vapor_core_cs_a) * Dz((m_dot_v(z) * m_dot_v(z) / (rho_v * vapor_core_cs_a)))

momentum_term = (2 * m_dot_v(z) / (vapor_core_cs_a^2 * rho_v)) * m_dot_per_length(z)

eq_P_v = Dz(P_v(z)) ~ - (friction_term) - momentum_term

friction_term_l = (mu_l / (wick_permeability * wick_cs_a * rho_l)) * m_dot_l(z)
gravity_term = rho_l * g * sin(pipe_tilt_angle)

eq_P_l = Dz(P_l(z)) ~ - (friction_term_l) + gravity_term

#Thermodynamic Equilibrium
eq_thermo = P_v(z) ~ P_ref * exp((H_vap_l / fluid_specific_gas_constant) * (1/T_ref - 1/T_v(z))) 

#END OF EQUATION DEFINITIONS 

z_start = 0
z_end = pipe_length
"""
max_capillary_P = 2 * fluid_surface_tension / wick_inner_radius 
delta_P_v = P_v(z=z_start) - P_v(z)
delta_P_l = P_l(z) - P_l(z_end)
grav_head = liquid_fluid_rho * g * pipe_length * sin(pipe_tilt_angle)
"""
eqs = [eq_energy_balance, eq_axial_conduction, eq_m_dot_v, eq_m_dot_l, eq_P_v, eq_P_l, eq_thermo]

domains = [z ∈ Interval(z_start, z_end)]
bcs = [
    m_dot_v(z_start) ~ 0,
    m_dot_v(z_end) ~ 0,
    #m_dot_v(z_end) ~ 0,
    m_dot_l(z_start) ~ 0,
    m_dot_v(z_end) ~ 0,
    #m_dot_l(z_end) ~ 0,
    # Ends are insulated (no axial heat flux)
    q_axial(z_start) ~ 0,
    q_axial(z_end) ~ 0,
    # Set a reference pressure. P_l=0 at the start of the liquid return path.
    P_l(z_start) ~ P_ref, # Set to a reference like 1 atm in Pa
    P_l(z_end) ~ P_ref, # Set to a reference like 1 atm in Pa
    P_v(z_start) ~ P_ref, # Set to a reference like 1 atm in Pa
    P_l(z_end) ~ P_ref, # Set to a reference like 1 atm in Pa
]

@named pdesystem = PDESystem(eqs, bcs, domains, [z], [T_v(z) P_v(z) P_l(z) m_dot_v(z) m_dot_l(z) q_axial(z) m_dot_per_length(z)])

z_points = 51
discretization = MOLFiniteDifference([z => z_points], nothing) #advection_scheme=(something) can either be Upwind or WENO
@time prob = discretize(pdesystem, discretization)

sol = solve(prob) #should probably implement callback = (something) here


