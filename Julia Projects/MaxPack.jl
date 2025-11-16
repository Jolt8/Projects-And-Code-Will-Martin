
""" Brainstorming
The Idea: 
    A julia toolkit for modelling the packing efficiency of different types of packing media 
    Honestly, it could also be a freecad or blender addon

Notes
    I could either use a game engine based physics solver or a physics based physics solver 

Features
    - Ability to import either an STL or mesh of the container and packing media
    - To improve performance, I could have it where the simulation stops adding packing media and then waits for a little and then freezes the simulation 
    and saves the state before adding more
    - At the end, once the maximum amount of packing is achieved, the state could become frozen and then exported to be used in a CFD program like openfoam
        - This would probably require a mesh of the packing media before the simulation is run
    - Uncertainty in the packing media's shape, 
        - For basic shapes, the scale and shape could vary a little 
        - For more complex shapes, the scale could vary or the user could provide multiple different types of files for the packing media and then randomly cycle between them
        - For the point above, this could also allow multiple types of packing media to be used by not randomly cycle in between them and choosing a predetermined distribution 
    - Ability to run seperate runs on different threads (maybe like 10 runs) to determine the most optimal packing

Outputs from a run
- Packing efficiency u"cm^3/cm^3"
- Final catalyst weight 
- catalyst surface area 


#Options for physics: Modia.jl (75 stars), RigidBodySim.jl (74 stars), Meshes.jl (436 stars), and writing my own solver (hard)

"""