# Author: Ivan Bioli (https://github.com/IvanBioli)

import Pkg
# Pkg.activate("elementifinitiunipv_pkg"; shared=false) # Uncomment this line if you want to activate the package environment
Pkg.add(name="Gridap", version="0.17")
Pkg.add(["GridapGmsh"])
Pkg.add(["PlotlyJS", "Plots", "LaTeXStrings"])
Pkg.add([
    "CairoMakie",
    "Gmsh",
    "Memoize",
    "Meshes",
    "Revise",
    "SparseArrays",
    "TriplotRecipes"
])