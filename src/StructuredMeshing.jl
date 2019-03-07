module StructuredMeshing

using LinearAlgebra

import Base.angle
import Base.show

export emptyMeshDef, addVertice, extrude, transitionextrude

include("definition.jl")
include("transfinite.jl")
include("generate.jl")
include("show.jl")

end