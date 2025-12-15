function calc_volume_integral!(du, u, mesh,
                               have_nonconservative_terms, equations,
                               volume_integral::VolumeIntegralWeakForm,
                               dg::DGSEM, cache)
    @threaded for element in eachelement(dg, cache)
        weak_form_kernel!(du, u, element, mesh,
                          have_nonconservative_terms, equations,
                          dg, cache)
    end

    return nothing
end

@inline function weak_form_kernel!(du, u,
                                   element, mesh::Union{TreeMesh{1}, StructuredMesh{1}},
                                   have_nonconservative_terms::False, equations,
                                   dg::DGSEM, cache, alpha = true)
    # true * [some floating point value] == [exactly the same floating point value]
    # This can (hopefully) be optimized away due to constant propagation.
    @unpack derivative_dhat = dg.basis

    for i in eachnode(dg)
        u_node = get_node_vars(u, equations, dg, i, element)

        flux1 = flux(u_node, 1, equations)
        for ii in eachnode(dg)
            multiply_add_to_node_vars!(du, alpha * derivative_dhat[ii, i], flux1,
                                       equations, dg, ii, element)
        end
    end

    return nothing
end

@inline function flux_differencing_kernel!(du, u,
                                           element,
                                           mesh::Union{TreeMesh{1}, StructuredMesh{1}},
                                           have_nonconservative_terms::False, equations,
                                           volume_flux, dg::DGSEM, cache, alpha = true)
    # true * [some floating point value] == [exactly the same floating point value]
    # This can (hopefully) be optimized away due to constant propagation.
    @unpack derivative_split = dg.basis

    # Calculate volume integral in one element
    for i in eachnode(dg)
        u_node = get_node_vars(u, equations, dg, i, element)

        # All diagonal entries of `derivative_split` are zero. Thus, we can skip
        # the computation of the diagonal terms. In addition, we use the symmetry
        # of the `volume_flux` to save half of the possible two-point flux
        # computations.

        # x direction
        for ii in (i + 1):nnodes(dg)
            u_node_ii = get_node_vars(u, equations, dg, ii, element)
            flux1 = volume_flux(u_node, u_node_ii, 1, equations)
            multiply_add_to_node_vars!(du, alpha * derivative_split[i, ii], flux1,
                                       equations, dg, i, element)
            multiply_add_to_node_vars!(du, alpha * derivative_split[ii, i], flux1,
                                       equations, dg, ii, element)
        end
    end
end

@inline function multiply_add_to_node_vars!(u, factor, u_node, equations, solver::DG,
                                            indices...)
    for v in eachvariable(equations)
        u[v, indices...] = u[v, indices...] + factor * u_node[v]
    end
    return nothing
end