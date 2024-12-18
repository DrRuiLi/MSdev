setClass("Molecule_igraph",
         slots = list(
           molecule_info = "list",
           sdf = "SDF",
           igraph = "ANY",
           isotopomer = "data.frame"
         ))
