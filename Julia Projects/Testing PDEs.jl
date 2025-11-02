using Unitful, DifferentialEquations
using Clapeyron
using Plots

using ModelingToolkit
using MethodOfLines
using DomainSets
# Parameters, variables, and derivatives
@parameters t x
@parameters dS dI brn ϵ
@variables S(..) I(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

# Define functions
function γ(x)
    y = x + 1.0
    return y
end

function ratio(x, brn, ϵ)
    y = brn + ϵ * sin(2 * pi * x)
    return y
end

# 1D PDE and boundary conditions
eq = [Dt(S(t, x)) ~ dS * Dxx(S(t, x)) - ratio(x, brn, ϵ) * γ(x) * S(t, x) * I(t, x) / (S(t, x) + I(t, x)) + γ(x) * I(t, x),
    Dt(I(t, x)) ~ dI * Dxx(I(t, x)) + ratio(x, brn, ϵ) * γ(x) * S(t, x) * I(t, x) / (S(t, x) + I(t, x)) - γ(x) * I(t, x)]
bcs = [S(0, x) ~ 0.9 + 0.1 * sin(2 * pi * x),
    I(0, x) ~ 0.1 + 0.1 * cos(2 * pi * x),
    Dx(S(t, 0)) ~ 0.0,
    Dx(S(t, 1)) ~ 0.0,
    Dx(I(t, 0)) ~ 0.0,
    Dx(I(t, 1)) ~ 0.0]

# Space and time domains
domains = [t ∈ Interval(0.0, 10.0),
    x ∈ Interval(0.0, 1.0)]

# PDE system
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [S(t, x), I(t, x)], [dS => 0.5, dI => 0.1, brn => 3, ϵ => 0.1])

# Method of lines discretization
# Need a small dx here for accuracy
dx = 0.01
order = 2
discretization = MOLFiniteDifference([x => dx], t)

# Convert the PDE problem into an ODE problem
@time prob = discretize(pdesys, discretization);