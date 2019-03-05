struct Boundary
    boundtype :: Symbol
    vertice_start :: Int64
    vertice_end :: Int64
    node_num :: Int64
    bias:: Float64
    
    function Boundary(boundtype::Symbol, vertice_start::Int64, vertice_end::Int64, node_num::Int64, bias::Float64=1.0)
        new(boundtype, vertice_start, vertice_end, node_num, bias)
    end
end

struct MeshDef
    blocks :: Vector{Dict{Symbol, Any}}
    bounds :: Vector{Boundary}
    vertices :: Vector{Vector{Float64}}
end

function emptyMeshDef()
    blocks = Vector{Dict{Symbol,Any}}()
    bounds = Vector{Boundary}()
    vertices = Vector{Vector{Float64}}()
    
    MeshDef(blocks, bounds, vertices)
end

function addVertice(meshdef::MeshDef, coordinates::Vector{Float64})
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

function extrude(meshdef::MeshDef, vert_start::Int64, direction::Vector{Float64}, length::Float64, node_num::Int64, bias::Float64=1.0)
    dir_norm = direction/norm(direction)
    
    coords_start = meshdef.vertices[vert_start]
    coords_new = coords_start + dir_norm * length
    
    vert_end = addVertice(meshdef, coords_new)
    
    bound_tmp = Boundary(:straight, vert_start, vert_end, node_num, -bias)
    
    push!(meshdef.bounds, bound_tmp)
    
    bound_id = size(meshdef.bounds, 1)
    
    return BoundaryLink(bound_id, vert_start, vert_end)
end

function connect(meshdef::MeshDef, vert_start::Int64, vert_end::Int64, node_num::Int64, bias::Float64=1.0)
    
    push!(meshdef.bounds, Boundary(:straight, vert_start, vert_end, node_num, bias))
    
    bound_id = size(meshdef.bounds, 1)
    
    return BoundaryLink(bound_id, vert_start, vert_end)
end

function extrude(meshdef::MeshDef, boundlink::BoundaryLink, direction::Vector{Float64}, length::Float64, node_num::Int64, bias::Float64=1.0;
        blocktype::Symbol=:transfinite, parallel_node_num::Int64=0)
    
    bounds = Vector{Int64}()
    
    boundlink1 = extrude(meshdef, boundlink.vert_start, direction, length, node_num, bias)
    
    boundlink2 = extrude(meshdef, boundlink.vert_end, direction, length, node_num, bias)
    
    startbound = meshdef.bounds[boundlink.id]
    
    # parallel_node_num is for transition mesh
    if parallel_node_num == 0
        parallel_node_num = startbound.node_num
    end
    
    boundlink3 = connect(meshdef, boundlink2.vert_end, boundlink1.vert_end, parallel_node_num, -startbound.bias)        
    
    bounds = [boundlink.id, boundlink2.id, boundlink3.id, -boundlink1.id]
    
    blockdef = Dict(:type=>blocktype, :bounds=>bounds)
    
    push!(meshdef.blocks, blockdef)
    
    block_id = size(meshdef.blocks, 1)
    
    BlockLink(block_id, boundlink3, boundlink1, boundlink2, 
        boundlink1.vert_start, boundlink1.vert_end)
end

function extrude(meshdef::MeshDef, boundlink::BoundaryLink, integrate::Vector{BoundaryLink})
    error("Not Implemented yet!")
end

function angle(vec_a::Vector{Float64}, vec_b::Vector{Float64})
    acos((vec_a'*vec_b) / (norm(vec_a) * norm(vec_b)))
end

function transitionextrude(meshdef::MeshDef, boundlink::BoundaryLink, direction::Vector{Float64}, layers::Int64)
    
    bound_start = meshdef.bounds[boundlink.id]
    
    if bound_start.boundtype != :straight && bound_start.bias != 1.0
        error("Transition block must start from straight, equidistant boundary!")
    end
    
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
    blocklink = extrude(meshdef, boundlink, vec, distance, layers+1, 0.5, blocktype=:transition,
                        parallel_node_num = N)
    
    return blocklink    
end


function defineCartesian(xmin::Float64, ymin::Float64, xmax::Float64, ymax::Float64, nx::Int64, ny::Int64)
    # Create vertices
    vertices = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]]
    
    # Create Bounds
    bounds = Vector{Boundary}(undef, 4)
    for i in 1:4
        n = i in (1,3) ? nx : ny
        bounds[i] = Boundary(:equidistant, i, i%4+1, n)
    end

    blocks = [Dict(:type=>:cartesian, :bounds=>[1,2,3,4])]
        
    return MeshDef(blocks, bounds, vertices)
end