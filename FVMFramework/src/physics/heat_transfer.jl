function get_k_effective(k_a, k_b)
    return 2 * k_a * k_b / (k_a + k_b)
end

function numerical_flux(k_avg, T_L, T_R, area, dist)
    grad_T = (T_R - T_L) / dist
    q = -k_avg * grad_T
    return q * area
end

function diffusion_temp_exchange!(
        #mutated vars
        du_temp_a, du_temp_b,
        #u data
        temp_a, temp_b,
        #geometry data
        connection_area, connection_distance,
        #props
        k_a, k_b,
        #other data

    )
    k_effective = get_k_effective(k_a, k_b)

    F = numerical_flux(k_effective, temp_a, temp_b, connection_area, connection_distance)
    
    du_temp_a -= F
    du_temp_b += F
end