"""
Stuctured Mesh Generation package. 
Provides tools for structured mesh definition
generation and export to abaqus.

Copyright 2019 Jonathan Hahn

Licensed under MIT License, see LICENSE.md
"""
module StructuredMeshing

using LinearAlgebra

import Base.angle
import Base.show

# types
export MeshDef, Mesh

# functions

## mesh definition
export emptyMeshDef, addVertice, extrude, transitionextrude, defineCartesian, connect

## mesh generation
export mesh

## abaqus export
export writeAbq

include("definition.jl")
include("transfinite.jl")
include("generate.jl")
include("show.jl")

end