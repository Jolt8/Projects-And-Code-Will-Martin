function species_numerical_flux(rho_avg, diffusion_coeff, mass_fraction_a, mass_fraction_b, area, dist) #not sure what k_average would represent in respect to diffusion
    concentration_gradient = (mass_fraction_b - mass_fraction_a) / dist
    diffusion = -rho_avg * diffusion_coeff * concentration_gradient
    return diffusion * area
end

function diffusion_mass_fraction_exchange!(
        #mutated vars
        du_species_mass_fractions_a, du_species_mass_fractions_b,
        #u data
        species_mass_fractions_a, species_mass_fractions_b,
        #geometry data
        area, dist,
        #props
        diffusion_coeff_a, diffusion_coeff_b,
        species_diffusion_coeffs_a, species_diffusion_coeffs_b
        #other data
    )

    rho_avg = 0.5 * (rho_a + rho_b)

    for i in eachindex(species_mass_fractions_a)
        diffusion_coeff_effective = harmonic_mean(species_diffusion_coeffs_a[i], species_diffusion_coeffs_b[i])
        numerical_flux = species_numerical_flux(rho_avg, diffusion_coeff_effective, species_mass_fractions_a[i], species_mass_fractions_b[i], area, dist)
        du_species_mass_fractions_a -= numerical_flux
        du_species_mass_fractions_b += numerical_flux
    end
end