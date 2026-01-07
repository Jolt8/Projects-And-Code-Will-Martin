using Unitful
using Symbolics, Nemo, Groebner
using Clapeyron
using NLsolve
using ModelingToolkit
# Process Variables (Temperatures, Flow, Heat)
@variables Q, LMTD, F_correction
@variables T_h_in, T_h_out, T_c_in, T_c_out
@variables m_h, m_c
@variables Cp_h, Cp_c

# Geometry Variables
@variables A, d_i, d_o, L_tube, N_tubes, N_passes, B_spacing, N_baffles, D_shell, P_t
@variables wall_thickness, k_wall

# Heat Transfer Coefficients & Resistances
@variables U, h_i, h_o
@variables Rf_i, Rf_o  # Fouling factors
@variables k_fluid_t, k_fluid_s
@variables mu_t, mu_s, mu_w_s  # Viscosity (tube, shell, wall)
@variables rho_t, rho_s

# Dimensionless Numbers & Intermediate Calcs
@variables Re_t, Pr_t, Nu_t, v_t
@variables Re_s, Pr_s, Nu_s, D_e, A_s, G_s
@variables n_dittus # Exponent for Dittus-Boelter (0.4 heating, 0.3 cooling)

# Pressure Drop Variables
@variables dP_t, dP_s
@variables f_t, f_s # Friction factors

# 1. Heat Duty (Energy Balance)
# ------------------------------------------
eq_energy_h = Q ~ m_h * Cp_h * (T_h_in - T_h_out)
eq_energy_c = Q ~ m_c * Cp_c * (T_c_out - T_c_in)

# 2. Rate Equation & LMTD
# ------------------------------------------
# Intermediate delta T variables for LMTD readability
@variables dt1, dt2

eq_dt1 = dt1 ~ T_h_in - T_c_out
eq_dt2 = dt2 ~ T_h_out - T_c_in

eq_LMTD = LMTD ~ (dt1 - dt2) / log(dt1 / dt2)

eq_rate = Q ~ U * A * F_correction * LMTD

# 3. Overall Heat Transfer Coefficient (Resistance Sum)
# ------------------------------------------
eq_U_inv = 1/U ~ (1/h_o) + Rf_o + 
                 (d_o * log(d_o/d_i)) / (2 * k_wall) + 
                 Rf_i * (d_o/d_i) + 
                 (1/h_i) * (d_o/d_i)

# 4. Tube-Side Heat Transfer (Dittus-Boelter Correlation)
# ------------------------------------------
# Area flow relationship for velocity
eq_v_t = v_t ~ m_h / (rho_t * (N_tubes/N_passes) * (pi * (d_i^2) / 4))

# Reynolds, Prandtl, Nusselt
eq_Re_t = Re_t ~ (rho_t * v_t * d_i) / mu_t
eq_Pr_t = Pr_t ~ (Cp_h * mu_t) / k_fluid_t
eq_Nu_t = Nu_t ~ 0.023 * (Re_t^0.8) * (Pr_t^n_dittus)

# Solve for h_i
eq_h_i = h_i ~ (Nu_t * k_fluid_t) / d_i

# 5. Shell-Side Heat Transfer (Kern Method)
# ------------------------------------------
# Equivalent Diameter (Assuming Square Pitch Layout here)
eq_De = D_e ~ (4 * ((P_t^2) - (pi * (d_o^2) / 4))) / (pi * d_o)

# Cross-flow area at shell equator
eq_As = A_s ~ ((P_t - d_o) * D_shell * B_spacing) / P_t

# Mass velocity and Reynolds
eq_Gs = G_s ~ m_c / A_s
eq_Re_s = Re_s ~ (D_e * G_s) / mu_s
eq_Pr_s = Pr_s ~ (Cp_c * mu_s) / k_fluid_s

# Nusselt (Kern)
eq_Nu_s = Nu_s ~ 0.36 * (Re_s^0.55) * (Pr_s^(1//3)) * ((mu_s / mu_w_s)^0.14)

# Solve for h_o
eq_h_o = h_o ~ (Nu_s * k_fluid_s) / D_e

# 6. Hydraulic / Pressure Drop
# ------------------------------------------
# Tube-side Pressure Drop (Friction + Return losses)
# Approximation: f_t = 0.046 * Re^-0.2 for smooth pipes (turbulent)
eq_f_t = f_t ~ 0.046 * (Re_t^(-0.2)) 
eq_dP_t = dP_t ~ ((4 * f_t * ((L_tube * N_passes) / d_i)) + (4 * N_passes)) * ((rho_t * (v_t^2)) / 2)

# Shell-side Pressure Drop
# Approximation: f_s = exp(0.576 - 0.19*ln(Re_s)) is common, 
# but here we keep f_s as a variable to be solved or supplied
eq_dP_s = dP_s ~ (f_s * (G_s^2) * D_shell * (N_baffles + 1)) / (2 * rho_s * D_e * ((mu_s / mu_w_s)^0.14))

# 7. Geometric Closure
# ------------------------------------------
eq_Area_geo = A ~ N_tubes * pi * d_o * L_tube
eq_Nb = N_baffles ~ (L_tube / B_spacing) - 1

# ==========================================
# SYSTEM ASSEMBLY
# ==========================================
# You can collect these into a system array for processing
all_eqs = [
    eq_energy_h,
    eq_energy_c,
    eq_dt1,
    eq_dt2,
    eq_LMTD,
    eq_rate,
    eq_U_inv,
    eq_v_t,
    eq_Re_t,
    eq_Pr_t,
    eq_Nu_t,
    eq_h_i,
    eq_De,
    eq_As,
    eq_Gs,
    eq_Re_s,
    eq_Pr_s,
    eq_Nu_s,
    eq_h_o,
    eq_f_t,
    eq_dP_t,
    eq_dP_s,
    eq_Area_geo,
    eq_Nb
]

try 
    symbolic_linear_solve(eq_h_o, h_o)
catch
    symbolic_solve(eq_h_o, h_o)[1]
end

#eq_Gs_substituted = substitute(eq_Gs, Dict(A_s => eq_As.rhs))

vars = [
    Q,              # Heat Duty
    LMTD,           # Log Mean Temp Difference
    T_h_out,        # Hot outlet temp
    T_c_out,        # Cold outlet temp
    dt1, dt2,       # LMTD intermediates
    U,              # Overall Heat Transfer Coeff
    h_i, h_o,       # Film coefficients
    v_t,            # Tube velocity
    Re_t, Pr_t, Nu_t, # Tube dimensionless numbers
    f_t,            # Tube friction factor
    dP_t,           # Tube pressure drop
    D_e,            # Equivalent diameter
    A_s,            # Shell cross-flow area
    G_s,            # Shell mass velocity
    Re_s, Pr_s, Nu_s, # Shell dimensionless numbers
    dP_s,           # Shell pressure drop
    A,              # Heat Transfer Area (Calculated)
    N_baffles       # Number of Baffles (Calculated)
]

# 2. Define the Parameters (Fixed Inputs)
# ---------------------------------------------------------
# These are the values you must provide to solve the system.
# Note: 'f_s' is listed here because there was no equation provided 
# to calculate it, so it must be a user input.
pars = [
    T_h_in, T_c_in,
    m_h, m_c,
    Cp_h, Cp_c,
    F_correction,
    A_s, d_i, d_o, L_tube, N_tubes, N_passes, B_spacing, D_shell, P_t,
    k_wall,
    Rf_i, Rf_o,
    k_fluid_t, k_fluid_s,
    mu_t, mu_s, mu_w_s,
    rho_t, rho_s,
    n_dittus,
    f_s             # Treated as parameter (input) as no correlation eq provided
]

# 3. Create the System
# ---------------------------------------------------------
# We pass the list of equations, the unknowns, and the parameters.
@named sys = NonlinearSystem(all_eqs, vars, pars)

# 4. Simplify the System (Recommended)
# ---------------------------------------------------------
# This performs the substitution you asked about in the previous step automatically.
# It reduces the 24 equations down to just the core nonlinear ones 
# (likely just the Energy Balance + Rate Equation loop).
simple_sys = structural_simplify(sys)

# You can now inspect what is left to solve:
println("Original equations: ", length(equations(sys)))
println("Simplified equations: ", length(equations(simple_sys)))