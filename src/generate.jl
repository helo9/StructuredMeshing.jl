mutable struct Mesh
    nodes :: Vector{Vector{Float64}}
    elements :: Vector{Vector{Int64}}
end

function emptyMesh()
    Mesh(Vector{Float64}[], Vector{Int64}[])
end

function getBiased(boundary::Boundary)
    N = boundary.node_num-1
    α = abs(boundary.bias)

    us = [0.0, [(α.^i-1)/(α^N-1) for i in range(0, stop=N, step=1)]...]

    if boundary.bias > 0
        return us
    else
        return 1 .- us
    end
end

function meshBoundary!(mesh::Mesh, boundary::Boundary, vertices, skipFirst::Bool, skipLast::Bool)
    if boundary.boundtype != :straight
        throw(DomainError())
    end
    
    # get coordinates of start and end point
    coords1 = vertices[boundary.vertice_start]
    coords2 = vertices[boundary.vertice_end]

    # Stretching functions
    Δcoords = coords2 - coords1
    N = boundary.node_num
    is = collect(range(0, stop=N-1, step=1))
    
    if abs(boundary.bias) == 1.0
        N = boundary.node_num
        nodes = (coords1 .+ Δcoords * (is ./ (N-1))')'
    else
        nodes = coords1' .+ Δcoords' .* getBiased(boundary)
    end
        
    #
    id1 = skipFirst ? 2 : 1
    id2 = skipLast ? 1 : 0
    
    node_ids = appendNodes!(mesh, nodes[id1:end-id2,:])
    
    return node_ids
end

function meshBoundaries!(mesh::Mesh, meshdefinition::MeshDef)
    # create Dict for Node mapping, mapping of node_ids to boundary respectiveley vertice ids
    nodeMapping = Dict(:vertice => Dict{Int64, Int64}(), :boundary => Dict{Int64, Vector{Int64}}())
    
    # create nodes on bounds
    for (boundary_id, boundary) in enumerate(meshdefinition.bounds)
        
        # when first or last node (at vertex) is already created with other boundary -> skip
        skipFirst = skipLast = false
        if verthasNode(nodeMapping, boundary.vertice_start)
            skipFirst = true
        end
        
        if verthasNode(nodeMapping, boundary.vertice_end)
            skipLast = true
        end
        
        # create Boundary Notes and append them to mesh
        node_ids = meshBoundary!(mesh, boundary, meshdefinition.vertices, skipFirst, skipLast)
        
        if sign(boundary_id) == -1
            node_ids = node_ids[end:-1:1]
        end
        
        nodeMapping[:boundary][abs(boundary_id)] = node_ids
        
        # add created nodes to nodeMapping
        if !skipFirst
            nodeMapping[:vertice][boundary.vertice_start] = node_ids[1]
        else
            insert!(nodeMapping[:boundary][abs(boundary_id)], 1, nodeMapping[:vertice][boundary.vertice_start])
        end
        if !skipLast
            nodeMapping[:vertice][boundary.vertice_end] = node_ids[end]
        else
            push!(nodeMapping[:boundary][abs(boundary_id)], nodeMapping[:vertice][boundary.vertice_end])
        end
    end
            
    return nodeMapping
end

function mesh(meshdefinition::MeshDef)
    
    # disable garbage collection during mesh creation
    GC.enable(false)
    
    # create an empty Mesh object
    mesh = emptyMesh()
    
    # generate nodes on boundaries
    nodeMapping = meshBoundaries!(mesh, meshdefinition)
    
    # generate blocks
    for (block_id, block) in enumerate(meshdefinition.blocks)                
        meshBlock!(mesh, block, meshdefinition.bounds, nodeMapping)
    end
    
    # enable garbage collection again
    GC.enable(true)
    
    return mesh
end

function verthasNode(nodeMapping, vert_id)
    return haskey(nodeMapping[:vertice], vert_id)
end
        
function getBoundaryNodeIds(bound_id, nodeMapping)
    boundarynodes = nodeMapping[:boundary][abs(bound_id)]
    
    # sign is used as direction
    if bound_id > 0
        return boundarynodes
    else
        return boundarynodes[end:-1:1]
    end
end

function meshBlock!(mesh::Mesh, block, bounds, nodeMapping)
    if block[:type] == :transfinite
        meshTransfiniteBlock!(mesh, block, bounds, nodeMapping)
    elseif block[:type] == :transition
        meshTransitionBlock!(mesh, block, bounds, nodeMapping)
    else
        blocktype_str = string(block[:type])
        throw(UndefVarError("Not defined blocktype ($blocktype_str)"))
    end
                    
end
                
function meshTransfiniteBlock!(mesh::Mesh, block, bounds, nodeMapping)
    
    # extract blocks boundaries
    boundary_ids = block[:bounds]
    boundaries = bounds[abs.(boundary_ids)]
        
    # number of nodes in directions 1 (bound1, -bound3)
    #   and 2 (bound2, -bound4)
    n1 = boundaries[1].node_num
    n2 = boundaries[2].node_num
    # TODO: check other bounds, throw Exception
    
    # get boundary nodes for first boundary from nodeMapping
    node_ids = getBoundaryNodeIds(boundary_ids[1], nodeMapping)
    
    # coords of second bound
    coords = mesh.nodes[getBoundaryNodeIds(boundary_ids[2], nodeMapping)]
                 
    # iterate over second dimension of block
    for i in 2:n2
        
        # store recent nodes
        last_node_ids = node_ids
        
        # generate next nodes
        if i<n2
            # calculate translation
            translation = coords[i]-coords[i-1]

            # copy nodes
            new_node_ids = translateNodes!(mesh, node_ids[2:end-1], translation, copynodes=true)
            
            node_ids = [getBoundaryNodeIds(boundary_ids[4], nodeMapping)[end-(i-1)],
                        new_node_ids..., 
                        getBoundaryNodeIds(boundary_ids[2], nodeMapping)[i]] 
        else
            node_ids = getBoundaryNodeIds(boundary_ids[3], nodeMapping)[end:-1:1]
        end
            
        for j in 2:n1
            # create Elements
                
            el_node1 = last_node_ids[j-1]
            el_node2 = last_node_ids[j]
            el_node3 = node_ids[j]
            el_node4 = node_ids[j-1]
            
            el_nodes = [el_node1, el_node2, el_node3, el_node4]
            
            el_id = appendElement!(mesh, el_nodes)
        end
    end
end

function meshTransfiniteBlock2!(mesh::Mesh, block, bounds, nodeMapping)
    
    # extract blocks boundaries
    boundary_ids = block[:bounds]
    boundaries = bounds[abs.(boundary_ids)]

    # extract vertice coordinates
    firstvert(boundary_id) = boundary_id>0 ? bounds[boundary_id].vertice_start : bounds[-boundary_id].vertice_end
    vertice_ids = [firstvert(boundary_id) for boundary_id in boundary_ids]
    vertices = Dict(id=>mesh.nodes[nodeMapping[:vertice][id]] for id in vertice_ids)

    # define corner points
    P12, P14, P34, P32 = vertices[vertice_ids]
    
    # extract edge definitions
    linefunc(p1, p2) = (u) -> p1 .+ (p2.-p1) .* u
    c1 = linefunc(P12, P14)
    c2 = linefunc(P32, P12)
    c3 = linefunc(P32, P34)
    c4 = linefunc(P14, P34)
    
    # calculate node positions
    us = (α.^i-1)/(α^N_-1)
    calculateTransfiniteNodes(c1, c2, c3, c4, P12, P14, P34, P32, us, vs)
end
                            
function meshTransitionBlock!(mesh::Mesh, block, bounds, nodeMapping)
    println("Does nothing at all!")
end
                        
function translateNodes!(mesh::Mesh, node_ids::Vector{Int64}, translation::Vector{Float64}; copynodes::Bool=true)
    if !copynodes
        throw(DomainError())
    end
    
    coords = Vector{Vector{Float64}}()
        
    for node_id in node_ids
        push!(coords, mesh.nodes[node_id]+translation)
    end
            
    node_ids_new = appendNodes!(mesh, coords)
    
    return node_ids_new
end
                
function appendNodes!(mesh::Mesh, nodecoords::Vector{Vector{Float64}})
    start_id = size(mesh.nodes, 1)+1
    length = size(nodecoords, 1)
    
    append!(mesh.nodes, nodecoords)
    
    return collect(start_id:1:start_id+length-1)
end

function appendNodes!(mesh::Mesh, nodecoords::Array{Float64,2})
    start_id = size(mesh.nodes, 1)+1
    length = size(nodecoords, 1)

    for i in range(1, stop=length, step=1)
        push!(mesh.nodes, nodecoords[i,:])
    end

    return collect(start_id:1:start_id+length-1)
end

function appendElement!(mesh::Mesh, node_ids)
    push!(mesh.elements, node_ids)
end
                
