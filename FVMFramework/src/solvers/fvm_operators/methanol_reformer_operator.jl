
function methanol_reformer_f!(
        du, u, p, t,
        cell_neighbor_map,
        cell_volumes, cell_centroids, 
        connection_areas, connection_normals, 
        #cell volumes and cell centroids are accessed at the id of the cell
        connection_distances, unconnected_areas,
        #connection areas, normals, and distances are simply accessed by their location in the 
        #list which corresponds to the respective connection in cell_neighbor_map
        species_molecular_weights,
        cell_props_id_map, bc_sys::BoundarySystem, chem_phys::Vector{ChemPhysics}, #heat_phys::Vector{HeatPhysics}, 
        molar_concentrations_cache, net_rates_cache,
        ax, n_reactions, n_species
    )

    #A_Ea_pairs = eachcol(reshape(p, :, n_reactions))
    #unflattened_p would be [[reaction_1_kf_A, reaction_1_kf_Ea], [reaction_2_kf_A, reaction_2_kf_Ea], etc..] 

    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)
    
    du .= 0.0

    #connections loop
    for (i, (idx_a, idx_b)) in enumerate(cell_neighbor_map)
        prop_a = cell_props_id_map[idx_a]
        prop_b = cell_props_id_map[idx_b]

        k_a = chem_phys[prop_a].k #maybe use heat_phys later
        k_b = chem_phys[prop_b].k

        connection_area = connection_areas[i]
        connection_distance = connection_distances[i]

        diffusion_temp_exchange!(
            du.temp[idx_a], du.temp[idx_b],
            u.temp[idx_a], u.temp[idx_b],
            connection_area, connection_distance,
            k_a, k_b,
        )
    end

    molar_concentrations_cache = zeros(eltype(u.mass_fractions), length(u.mass_fractions[:, 1])) #just using mass fractions for cell 1, this may cause some issues later!
    net_rates_cache = zeros(eltype(u.mass_fractions), length(chem_phys[1].chemical_reactions))

    # Source and Capacity Loop
    for cell_id in bc_sys.free_idxs
        #basic props
        vol = cell_volumes[cell_id]
        props = cell_props_id_map[cell_id]

        k = chem_phys[props].k
        rho = chem_phys[props].rho
        cp  = chem_phys[props].cp
        cell_chemical_reactions_vec = chem_phys[props].chemical_reactions



        #chemical reactions loop
        #species_mass_fractions = u.mass_fractions[:, cell_id]
        cell_temp = u.temp[cell_id]
        species_mass_fractions = view(u.mass_fractions, :, cell_id) #we should maybe use views here, probably does matter that much
        
        react_cell!(
            @view(du.mass_fractions[:, cell_id]), @view(du.temp[cell_id:cell_id]), 
            molar_concentrations_cache, net_rates_cache, 
            species_mass_fractions, cell_temp,
            vol,
            rho, 
            species_molecular_weights, cell_chemical_reactions_vec
        )

        # heat source and capacity loop
        S = chem_phys[props].heat_vol_source_term * vol 
        # we should probably create separate containers in heat_phys for both source terms on a per area and per cell basis

        du.temp[cell_id] += S
        
        cap = rho * cp * vol
        du.temp[cell_id] /= cap
    end

    for cell_id in bc_sys.dirichlet_idxs
        du.temp[cell_id] = 0.0
    end
end