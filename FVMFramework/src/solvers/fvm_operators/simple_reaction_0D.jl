struct SimpleReactionBoundarySystem <: AbstractBoundarySystem
    boundary_map::MultiPhysicsBCs
    free_idxs::Vector{Int}
    
    temp_fixed_idxs::Vector{Int}
    chem_fixed_idxs::Vector{Int}
end
struct SimpleReactionPhysicsBCs <: MultiPhysicsBCs #this also defines the order of each variable in u used in the future
    temp_bcs::Vector{HeatBC}
    chem_bcs::Vector{ChemBC}
end

function simple_reaction_0D_f!(
        du, u, p, t,
        cell_neighbor_map,
        cell_volumes, cell_centroids, 
        connection_areas, connection_normals, connection_distances,
        #cell volumes and cell centroids are accessed at the id of the cell
        unconnected_areas,
        #connection areas, normals, and distances are simply accessed by their location in the 
        #list which corresponds to the respective connection in cell_neighbor_map
        species_molecular_weights,
        cell_props_id_map, bc_sys::AbstractBoundarySystem, chem_phys::Vector{SimpleChemPhysics}, #heat_phys::Vector{HeatPhysics}, 
        ax, n_reactions, n_species
    )

    u = ComponentVector(u, ax)
    du = ComponentVector(du, ax)
    
    du .= 0.0

    molar_concentrations_cache = zeros(eltype(u.mass_fractions), length(u.mass_fractions[:, 1])) #just using mass fractions for cell 1, this may cause some issues later!
    net_rates_cache = zeros(eltype(u.mass_fractions), length(chem_phys[1].chemical_reactions))

    # Source and Capacity Loop
    for cell_id in bc_sys.free_idxs
        #basic props
        vol = cell_volumes[cell_id]
        props = cell_props_id_map[cell_id]
        
        rho = chem_phys[props].rho
        cp  = chem_phys[props].cp
        reactions = chem_phys[props].chemical_reactions
        cell_kg_cat_per_m3_for_each_reaction = chem_phys[props].cell_kg_cat_per_m3_for_each_reaction


        #chemical reactions loop
        #mass_fractions = u.mass_fractions[:, cell_id]
        cell_temp = u.temp[cell_id]
        mass_fractions = view(u.mass_fractions, :, cell_id) #we should maybe use views here, probably does matter that much
        
        #react_cell! also adds the heat of reaction to the cell 
        react_cell!(
            @view(du.mass_fractions[:, cell_id]), @view(du.temp[cell_id:cell_id]), #fixed temp for now
            molar_concentrations_cache, net_rates_cache, 
            mass_fractions, cell_temp,
            vol,
            rho, 
            species_molecular_weights, reactions, cell_kg_cat_per_m3_for_each_reaction
        )

        # heat source and capacity loop
        S = chem_phys[props].heat_vol_source_term * vol 
        
        du.temp[cell_id] += S

        # ----- CAPACITY LOOP -----
        cell_mass = rho * vol
        du.temp[cell_id] /= (rho * vol * cp)
    end

    for cell_id in bc_sys.temp_fixed_idxs
        du.temp[cell_id] = 0.0
    end
end