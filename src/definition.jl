"""
Container for boundary definitions.
"""

abstract type AbstractBoundary end

struct StraightBoundary <: AbstractBoundary
    vertice_start :: Int64
    vertice_end :: Int64
    node_num :: Int64
    bias :: Float64
    
    function StraightBoundary(vertice_start::Int64, vertice_end::Int64, node_num::Int64, bias::Float64=1.0)
        new(vertice_start, vertice_end, node_num, bias)
    end
end
    
struct CircularBoundary <: AbstractBoundary
    vertice_start :: Int64
    vertice_end :: Int64
    node_num :: Int64
    bias :: Float64
    radius :: Float64
    
    function CircularBoundary(vertice_start::Int64, vertice_end::Int64, node_num::Int64, radius::Float64, bias::Float64=1.0)
        new(vertice_start, vertice_end, node_num, bias, radius)
    end
end

"""
Container for mesh defintion, consisting of multiple blocks.
"""
struct MeshDef
    blocks :: Vector{Dict{Symbol, Any}}
    bounds :: Vector{AbstractBoundary}
    vertices :: Vector{Vector{Float64}}
    default_elem_type :: String
end


function emptyMeshDef(default_elem_type="CPS4")
    blocks = Vector{Dict{Symbol,Any}}()
    bounds = Vector{AbstractBoundary}()
    vertices = Vector{Vector{Float64}}()
    
    MeshDef(blocks, bounds, vertices, "CPS4")
end

"""
    addVertex(meshdef, coordinates)

Adds a vertex to a mesh definition, returns the vertex's id
"""
function addVertex!(meshdef::MeshDef, coordinates::Vector{Float64})
    push!(meshdef.vertices, coordinates)
    
    return size(meshdef.vertices, 1)
end

struct BoundaryLink
    id :: Int64
    vert_start :: Int64
    vert_end :: Int64
end

struct BlockLink
    id :: Int64
    frontbound :: BoundaryLink
    leftbound :: BoundaryLink
    rightbound :: BoundaryLink
    vert_left :: Int64
    vert_right :: Int64
end

"""
    extrude(meshdef, vert_start, direction, length, node_num, bias)

Extrude an vertex - results in an boundary
"""
function extrude!(meshdef::MeshDef, vert_start::Int64, direction::Vector{Float64}, length::Float64, node_num::Int64, bias::Float64=1.0)
    dir_norm = direction/norm(direction)
    
    coords_start = meshdef.vertices[vert_start]
    coords_new = coords_start + dir_norm * length
    
    vert_end = addVertex!(meshdef, coords_new)
    
    bound_tmp = StraightBoundary(vert_start, vert_end, node_num, -bias)
    
    push!(meshdef.bounds, bound_tmp)
    
    bound_id = size(meshdef.bounds, 1)
    
    return BoundaryLink(bound_id, vert_start, vert_end)
end

"""
    connect(meshdef, vert_start, vert_end, node_num, bias)

Connect two vertices to get a boundary
"""
function connect!(meshdef::MeshDef, vert_start::Int64, vert_end::Int64, node_num::Int64, bias::Float64=1.0)
    
    push!(meshdef.bounds, StraightBoundary(vert_start, vert_end, node_num, bias))
    
    bound_id = size(meshdef.bounds, 1)
    
    return BoundaryLink(bound_id, vert_start, vert_end)
end

"""
    extrude(meshdef, boundlink, direction, length, node_num, bias)

Extrude boundary to build new block
"""
function extrude!(meshdef::MeshDef, boundlink::BoundaryLink, direction::Vector{Float64}, length::Float64, node_num::Int64, bias::Float64=1.0;
        blocktype::Symbol=:transfinite, parallel_node_num::Int64=0, element_type::String="")
    
    bounds = Vector{Int64}()
    
    boundlink1 = extrude!(meshdef, boundlink.vert_start, direction, length, node_num, bias)
    
    boundlink2 = extrude!(meshdef, boundlink.vert_end, direction, length, node_num, bias)
    
    startbound = meshdef.bounds[boundlink.id]
    
    # parallel_node_num is for transition mesh
    if parallel_node_num == 0
        parallel_node_num = startbound.node_num
    end
    
    boundlink3 = connect!(meshdef, boundlink2.vert_end, boundlink1.vert_end, parallel_node_num, -startbound.bias)        
    
    bounds = [boundlink.id, boundlink2.id, boundlink3.id, -boundlink1.id]
    
    if isempty(element_type)
        element_type = meshdef.default_elem_type
    end
    
    blockdef = Dict(:type=>blocktype, :bounds=>bounds, :elem_type=>element_type)
    
    push!(meshdef.blocks, blockdef)
    
    block_id = size(meshdef.blocks, 1)
    
    BlockLink(block_id, boundlink3, boundlink1, boundlink2, 
        boundlink1.vert_start, boundlink1.vert_end)
end

function extrude!(meshdef::MeshDef, boundlink::BoundaryLink, integrate::Vector{BoundaryLink})
    error("Not Implemented yet!")
end

function angle(vec_a::Vector{Float64}, vec_b::Vector{Float64})
    acos((vec_a'*vec_b) / (norm(vec_a) * norm(vec_b)))
end

"""
    transitionextrude(meshdef, boundlink, direction, layers)

Generate Transition Mesh from extrusion of boundary 
"""
function transitionextrude!(meshdef::MeshDef, boundlink::BoundaryLink, direction::Vector{Float64}, layers::Int64; element_type::String="")
    
    bound_start = meshdef.bounds[boundlink.id]
    
    # calculate Element Size on start boundary
    coords1 = meshdef.vertices[boundlink.vert_start]
    coords2 = meshdef.vertices[boundlink.vert_end]
    
    Δl = norm(coords2-coords1)/(bound_start.node_num-1)
    
    # calculate extrusion length (geometric series)
    #n = layers-1
    distance = 2 * Δl * (2^layers - 1)
    
    # calculate direction normal to start boundary
    normalvec = (coords2-coords1)/norm(coords2-coords1)
    
    vec = [-normalvec[2], normalvec[1]]
    
    if abs(angle(vec, direction)) >= pi/2
        vec *= -1
    end
    
    # reduced mesh density on created parallel boundary
    N = convert(Int64, floor((bound_start.node_num-1) / 2^layers))+1
    
    # do extrusion
    blocklink = extrude!(meshdef, boundlink, vec, distance, layers+1, -0.5, blocktype=:transition,
                        parallel_node_num = N, element_type=element_type)
    
    return blocklink    
end

"""
    rotate!(meshdef, vert_star, center, angle, node_num)

Generate new boundary by the rotation of a vertice around center by angle.
"""
function rotate!(meshdef, vert_start::Integer, center, angle::Real, node_num::Integer)
    coords_start = meshdef.vertices[vert_start]
    
    leg = coords_start - center
    
    radius = norm(leg)
    
    # rotate leg
    T = [cos(angle) -sin(angle); sin(angle) cos(angle)]
    leg_rot = T * leg
    
    coords_end = center + leg_rot
    
    vert_end = addVertex!(meshdef, coords_end)
    
    boundary = CircularBoundary(vert_start, vert_end, node_num, radius)
    
    push!(meshdef.bounds, boundary)
    
    bound_id = size(meshdef.bounds, 1)
    
    return BoundaryLink(bound_id, vert_start, vert_end)
end


"""
    rotate(meshdef, boundary, center, angle, fill, node_num)

Generate new block or boundary by rotation of block around center by angle.
"""
function rotate!(meshdef::MeshDef, boundary::BoundaryLink, center, angle::Real, node_num::Integer; fill=false, element_type::String="")
    boundlink1 = rotate!(meshdef, boundary.vert_start, center, angle, node_num)
    boundlink2 = rotate!(meshdef, boundary.vert_end, center, angle, node_num)
    
    boundlink3 = rotate_!(meshdef, meshdef.bounds[boundary.id], boundlink1.vert_end, boundlink2.vert_end, center, angle)
    
    bounds = [boundlink1.id, boundlink3.id, -boundlink2.id, -boundary.id]
    
    if isempty(element_type)
        element_tpye = meshdef.default_elem_type
    end
    
    blockdef = Dict(:type=>:transfinite, :bounds=>bounds, :elem_type=>meshdef.default_elem_type)
    
    push!(meshdef.blocks, blockdef)
    
    block_id = size(meshdef.blocks, 1)
    
    BlockLink(block_id, boundlink3, boundlink1, boundlink2, 
        boundlink1.vert_start, boundlink1.vert_end)
    
end

function rotate_!(meshdef::MeshDef, boundary::StraightBoundary, vert_start, vert_end, center, anlge::Real)
    
    boundary = StraightBoundary(vert_start, vert_end, boundary.node_num)
    
    push!(meshdef.bounds, boundary)
    
    bound_id = size(meshdef.bounds, 1)
    
    return BoundaryLink(bound_id, vert_start, vert_end)
end

"""
    defineCartesian(xmin, ymin, xmax, ymax, nx, ny)

Build up simple equidistant cartesian mesh
"""
function defineCartesian(xmin::Float64, ymin::Float64, xmax::Float64, ymax::Float64, nx::Int64, ny::Int64)
    # Create vertices
    vertices = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]]
    
    # Create Bounds
    bounds = Vector{Boundary}(undef, 4)
    for i in 1:4
        n = i in (1,3) ? nx : ny
        bounds[i] = StraightBoundary(i, i%4+1, n)
    end

    blocks = [Dict(:type=>:cartesian, :bounds=>[1,2,3,4], :elem_type=>"CPS4")]
        
    return MeshDef(blocks, bounds, vertices)
end