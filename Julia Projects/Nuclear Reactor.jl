using ModelingToolkit
using DifferentialEquations
using MethodOfLines
using DomainSets

@parameters v_1=2.0e7 v_2=2.2e5 D_1=1.5 D_2=0.4 Σa_1=0.01
@parameters Σa_2=0.08 νΣf_1=0.005 νΣf_2=0.09 Σs_12=0.02
@parameters β=0.0065 λ=0.075

@parameters t, r
@syms ϕ_1(t, r) ϕ_2(t, r) C(t, r)
Dt = Differential(t)
Dr = Differential(r)
Drr = Differential(r)^2

r_min = 0.0
r_max = 100.0
t_min = 0.0
t_max = 1.0

N = 100
dr = (r_min - r_max) / N

domains = [
    r ∈ Interval(r_min, r_max),
    t ∈ Interval(t_min, t_max)
]

fission_source = (1 - β) * (νΣf_1 * ϕ_1(t, r) + νΣf_2 * ϕ_2(t, r))

# Leakage term in 1D Cylindrical coordinates is (1/r) * D(r * D * ϕ)
# We can expand this using the product rule to D(D(ϕ)) + (1/r)*D(ϕ)
# This is often more numerically stable for MethodOfLines at r=0
leakage_1 = D_1 * (Drr(ϕ_1(t, r)) + (1/r) * Dr(ϕ_1(t, r)))
leakage_2 = D_2 * (Drr(ϕ_2(t, r)) + (1/r) * Dr(ϕ_2(t, r)))

# Removal term = Absorption + Scattering out of the group
removal_1 = (Σa_1 + Σs_12) * ϕ_1(t, r)
removal_2 = Σa_2 * ϕ_2(t, r)

# Source for group 2 is neutrons scattering down from group 1
source_2 = Σs_12 * ϕ_1(t, r)

eqs = [
    # Equation for fast neutrons (group 1)
    (1/v_1) * Dt(ϕ_1(t, r)) ~ leakage_1 - removal_1 + fission_source,

    # Equation for thermal neutrons (group 2)
    (1/v_2) * Dt(ϕ_2(t, r)) ~ leakage_2 - removal_2 + source_2 + λ * C(t, r),

    # Equation for delayed neutron precursors
    Dt(C(t, r)) ~ β * (νΣf_1 * ϕ_1(t, r) + νΣf_2 * ϕ_2(t, r)) - λ * C(t, r)
]


# ===================================================================
# 3. DEFINE BOUNDARY AND INITIAL CONDITIONS
# ===================================================================

# Initial Conditions: what the reactor looks like at t=0
# Let's start with a small, flat flux everywhere to kick off the reaction
initial_flux = 1.0e5
ic = [
    ϕ_1(0, r) ~ initial_flux,
    ϕ_2(0, r) ~ initial_flux,
    C(0, r) ~ (β / λ) * (νΣf_1 * initial_flux + νΣf_2 * initial_flux) # Assume equilibrium at t=0
]

# Boundary Conditions: what happens at the edges r=0 and r=r_max
bcs = [
    # At the center (r=0), zero gradient (symmetry)
    Dr(ϕ_1(t, 0)) ~ 0,
    Dr(ϕ_2(t, 0)) ~ 0,
    Dr(C(t, 0)) ~ 0,

    # At the outer edge (r=r_max), zero flux (vacuum)
    ϕ_1(t, r_max) ~ 0,
    ϕ_2(t, r_max) ~ 0,
    C(t, r_max) ~ 0
]

# Combine BCs and ICs
bcs = vcat(bcs, ic)

# ===================================================================
# 4. SOLVE THE SYSTEM
# ===================================================================

# Create the PDESystem
@named pdesys = PDESystem(eqs, bcs, domains, [t, r], [ϕ_1(t, r), ϕ_2(t, r), C(t, r)])

# Discretize the spatial domain (r) into N points
N = 101 # Number of grid points
discretization = MOLFiniteDifference([r => N], t)

# Convert the PDESystem into an ODEProblem
prob = discretize(pdesys, discretization)

# Solve the ODEProblem
# This is a stiff system, so we use a solver designed for them
println("Solving the system...")
sol = solve(prob, Rosenbrock23(), saveat=0.1)
println("Solution complete.")

# ===================================================================
# 5. VISUALIZE THE RESULTS
# ===================================================================

# Extract the independent and dependent variables from the solution
discrete_t = sol.t
discrete_r = sol.ivs[r]
sol_phi_1 = sol[ϕ_1(t, r)]
sol_phi_2 = sol[ϕ_2(t, r)]

# Plot the thermal flux (group 2) at different time points
p = plot(discrete_r, sol_phi_2[1, :], label="t = $(discrete_t[1]) s",
         xlabel="Radius (cm)", ylabel="Thermal Flux", title="Thermal Neutron Flux Profile")
for i in 2:length(discrete_t)
    plot!(p, discrete_r, sol_phi_2[i, :], label="t = $(discrete_t[i]) s")
end
display(p)