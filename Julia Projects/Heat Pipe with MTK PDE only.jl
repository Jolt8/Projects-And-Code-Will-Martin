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

#mu_l_func(T) = 2.415e-5 * 10^(248.6 / (T - 140.0)) # Liquid viscosity (Pa*s)
#mu_v_func(T) = 1.002e-6 + 1.986e-8*T + 4.58e-12*T^2 # Vapor viscosity (Pa*s)

#k_l_func(T) = 124.4 - 0.1138*T + 5.523e-5*T^2 - 1.184e-8*T^3 # Liquid conductivity (W/m-K)

#dynam_visc_l_func(T) = 0.2323 - 1.01e-4 * T 

# Wick and Wall properties
wall_k = ustrip(uconvert(u"W/(m*K)", 20.0u"W/(m*K)"))
wall_rho = ustrip(uconvert(u"kg/m^3", 8220u"kg/m^3"))
wall_Cp = ustrip(uconvertu"J/(kg*K)", 435u"J/(kg*K)")

saturated_wick_effective_k = ustrip(uconvert(u"W/(m*K)", 40u"W/(m*K)"))
wick_permeability = ustrip(uconvert(u"m^2", 1e-11u"m^2"))
wick_porosity = 0.6
g = ustrip(uconvert(u"m/s^2", 9.80665u"m/s^2"))
wick_effective_capillary_radius = ustrip(uconvert(u"m", 50u"μm")) #length usually in mm or μm

P_ref = ustrip(uconvert(u"Pa", 101325u"Pa"))
T_ref = ustrip(uconvert(u"K", 1156.0u"K"))
fluid_specific_gas_constant = ustrip(uconvert(u"J/(kg*K)", 361.3u"J/(kg*K)"))


@independent_variables t x
@variables T_w(..) T_l(..) P_l(..) u_l(..) T_v(..) P_v(..) rho_v(..) u_v(..) m_dot_evap_flux(..)
Dt = Differential(t)
Dx = Differential(x)


Q_total_in = 1000.0 # Total power input in Watts
q_evap = Q_total_in / evap_length
q_cond = -Q_total_in / cond_length
q_input(x) = ifelse(x <= evap_length, q_evap,
ifelse(x <= evap_length + adia_length, 0.0, q_cond))
@register_symbolic q_input(x)

h_wl = ustrip(uconvert(u"W/(m^2*K)", 1.0u"W/(m^2*K)"))
P_wl = 2 * pi * wick_outer_radius

wall_energy_eq = wall_rho * wall_Cp * wall_cs_a * Dt(T_w(t, x)) ~ Dx(wall_k * wall_cs_a * Dx(T_w(t, x))) + q_input(x)*wall_cs_a - h_wl * P_wl * (T_w(t, x) - T_l(t, x))

rho_l(t, x) = 1011.6 - 0.2208*T_l(t, x) - 1.46e-5*T_l(t, x)^2 + 5.69e-9*T_l(t, x)^3
Cp_l = ustrip(uconvert(u"J/(kg*K)", 1.25u"J/(g*K)")) #good estimate
H_vap_l(t, x) = (4.9886e6 - 5.751e2*T_l(t, x) - 2.158*T_l(t, x)^2) / 1000 # Latent heat (J/kg), converted from J/kmol
P_lv = 2 * pi * wick_inner_radius
dynam_visc_l(t, x) = 0.2323 - 1.01e-4 * T_l(t, x)

wick_energy_eq = (rho_l(t, x) * Cp_l * wick_porosity) * wick_cs_a * Dt(T_l(t, x)) + (rho_l(t, x) * Cp_l) * wick_cs_a * u_l(t, x) * Dx(T_l(t, x)) ~ Dx(saturated_wick_effective_k * wick_cs_a * Dx(T_l(t, x))) + h_wl * P_wl * (T_w(t, x) - T_l(t, x)) - m_dot_evap_flux(t, x) * P_lv * H_vap_l(t, x)

wick_mass_cons_eq = wick_porosity * wick_cs_a * Dx(u_l(t, x)) = - (m_dot_evap_flux(t, x) * P_lv) / rho_l

wick_momentum_eq = u_l(t, x) ~ - (wick_permeability / dynam_visc_l) * (Dx(P_l(t, x)) + rho_l * g * sin(pipe_tilt_angle))

rho_v(t, x) = P_v(t, x) / (361.3 * T_v(t, x))
vapor_mass_cons_eq = vapor_core_cs_a * Dt(rho_v(t, x)) + Dx(rho_v(t, x) * u_l(t, x) * vapor_core_cs_a) ~ m_dot_evap_flux(t, x) * P_lv

darcy_fric_factor = 0.01
vapor_momentum_eq = Dt(rho_v(t, x) * u_l(t, x)) + Dt(rho_v(t, x) * u_l(t, x)^2) ~ - Dx(P_v) - (darcy_fric_factor / (2 * 2 * sqrt(vapor_core_cs_a / pi))) * rho_v(t, x) 