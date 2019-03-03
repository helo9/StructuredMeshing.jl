mutable struct Mesh
    nodes :: Vector{Vector{Float64}}
    elements :: Vector{Vector{Int64}}
end

function emptyMesh()
    Mesh(Vector{Float64}[], Vector{Int64}[])
end

function appendBoundaryNodes!(mesh::Mesh, boundary::Boundary, vertices, skipFirst::Bool, skipLast::Bool)
    if boundary.boundtype != :straight
        throw(DomainError())
    end
    
    
    # get coordinates of start and end point
    coords1 = vertices[boundary.vertice_start]
    coords2 = vertices[boundary.vertice_end]

    # Stretching functions
    Δcoords = (coords2-coords1)
    N = boundary.node_num
    is = collect(range(0, stop=N-1, step=1))
    
    if abs(boundary.bias) == 1.0
        N = boundary.node_num
        Δcoords = coords2-coords1
        nodes = [coords1 + Δcoords * i/(N-1) for i in is]
    else
        α = abs(boundary.bias)
        N = boundary.node_num-1
        Δcoords *= sign(boundary.bias)
        
        startcoord = boundary.bias>0 ? coords1 : coords2
        nodes = [startcoord + Δcoords * (α.^i-1)/(α^N-1) for i in is]

    end
    
    #
    id1 = skipFirst ? 2 : 1
    id2 = skipLast ? 1 : 0
    
    nodeids = appendNodes!(mesh, nodes[id1:end-id2])
    
    return nodeids
end

function generateBoundaryMesh!(mesh::Mesh, meshdefinition::MeshDef)
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
        node_ids = appendBoundaryNodes!(mesh, boundary, meshdefinition.vertices, skipFirst, skipLast)
        
        nodeMapping[:boundary][boundary_id] = node_ids
        
        # add created nodes to nodeMapping
        if !skipFirst
            nodeMapping[:vertice][boundary.vertice_start] = node_ids[1]
        else
            insert!(nodeMapping[:boundary][boundary_id], 1, nodeMapping[:vertice][boundary.vertice_start])
        end
        if !skipLast
            nodeMapping[:vertice][boundary.vertice_end] = node_ids[end]
        else
            push!(nodeMapping[:boundary][boundary_id], nodeMapping[:vertice][boundary.vertice_end])
        end
    end
            
    return nodeMapping
end

function generateMesh(meshdefinition::MeshDef)
    
    # disable garbage collection during mesh creation
    GC.enable(false)
    
    # create an empty Mesh object
    mesh = emptyMesh()
    
    
    nodeMapping = generateBoundaryMesh!(mesh, meshdefinition)
    
    
    # generate blocks
    for (block_id, block) in enumerate(meshdefinition.blocks)
                
        block_bounds = [meshdef.bounds[abs(i)] for i in block[:bounds]]
                
        addBlock!(mesh, block, block_bounds, nodeMapping)
    end
    
    GC.enable(true)
            
    return mesh
end

function verthasNode(nodeMapping, vert_id)
    return haskey(nodeMapping[:vertice], vert_id)
end
        
function getBoundaryNodeIds(bound_id, nodeMapping)
    boundarynodes = nodeMapping[:boundary][abs(bound_id)]
    
    if bound_id > 0
        return boundarynodes
    else
        return boundarynodes[end:-1:1]
    end
end

function addBlock!(mesh::Mesh, block, bounds, nodeMapping)
    if !(block[:type] == :cartesian)
        throw(DomainError())
    end
    
    # extract blocks boundaries
    boundary_ids = block[:bounds]
    boundaries = bounds[boundary_ids]
    
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
            node_ids = nodeMapping[:boundary][boundary_ids[3]][end:-1:1]
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

function translateNodes!(mesh::Mesh, node_ids::Vector{Int64}, translation::Vector{Float64}; copynodes::Bool=true)
    if !copynodes
        throw(DomainError())
    end
    
    coords = []
        
    for node_id in node_ids
        push!(coords, mesh.nodes[node_id]+translation)
    end
            
    node_ids_new = appendNodes!(mesh, coords)
    
    return node_ids_new
end
                
function appendNodes!(mesh::Mesh, nodecoords)
    start_id = size(mesh.nodes, 1)+1
    length = size(nodecoords, 1)
    
    append!(mesh.nodes, nodecoords)
    
    return collect(start_id:1:start_id+length-1)
end

function appendElement!(mesh::Mesh, node_ids)
    push!(mesh.elements, node_ids)
end
                