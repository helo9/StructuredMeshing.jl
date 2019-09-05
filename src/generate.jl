import Statistics: mean

"""
Container for generated Mesh
"""
mutable struct Mesh
    nodes :: Vector{Vector{Float64}}
    elements :: Vector{Vector{Int64}}
end

function emptyMesh()
    Mesh(Vector{Float64}[], Vector{Int64}[])
end

function getFac(boundary::Boundary; flipdirection::Bool=false)
    
    N = boundary.node_num-1

    if abs(boundary.bias) != 1.0
        α = abs(boundary.bias)

        us = [[(α.^i-1)/(α^N-1) for i in range(0, stop=N, step=1)]...]

    else
        us = collect(0:1:N)/(N)
    end

    if (boundary.bias > 0) != (!flipdirection)
        return us
    else
        return (1 .- us)[end:-1:1]
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
    
    # distribution
    us = getFac(boundary)
    nodes = coords1' .+ Δcoords' .* us
        
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

"""
    mesh(meshdef)

Generate mesh from meshdefinition
"""
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
        meshTransfiniteBlock2!(mesh, block, bounds, nodeMapping)
    elseif block[:type] == :transition
        meshTransitionBlock!(mesh, block, bounds, nodeMapping)
    else
        blocktype_str = string(block[:type])
        throw(UndefVarError("Not defined blocktype ($blocktype_str)"))
    end
                    
end

function blockdef2fun(mesh, block, bounds, nodeMapping)
    # extract blocks boundaries
    boundary_ids = block[:bounds]
    boundaries = bounds[abs.(boundary_ids)]

    # extract vertice coordinates
    firstvert(boundary_id) = boundary_id>0 ? bounds[boundary_id].vertice_start : bounds[-boundary_id].vertice_end
    vertice_ids = [firstvert(boundary_id) for boundary_id in boundary_ids]
    vertices = Dict(id=>mesh.nodes[nodeMapping[:vertice][id]] for id in vertice_ids)

    println()
    println(vertice_ids)
    println(vertices)
    println(nodeMapping[:vertice])
    println()

    # define corner points
    P12, P14, P34, P32 = (vertices[i] for i in vertice_ids)
    
    # extract edge definitions
    linefunc(p1, p2) = (u) -> p1 .+ (p2.-p1) .* u
    c1 = linefunc(P12, P14)
    c2 = linefunc(P12, P32)
    c3 = linefunc(P32, P34)
    c4 = linefunc(P14, P34)

    return (P12, P14, P34, P32), (c1, c2, c3, c4)
end

function meshTransfiniteBlock2!(mesh::Mesh, block, bounds, nodeMapping)
    
    function addElements!(mesh, last_node_ids, node_ids)
        for j in 2:size(node_ids,1)
            el_node1 = last_node_ids[j-1]
            el_node2 = node_ids[j-1]
            el_node3 = node_ids[j]
            el_node4 = last_node_ids[j]

            el_nodes = [el_node1, el_node2, el_node3, el_node4]
            
            el_id = appendElement!(mesh, el_nodes)
        end
    end

    # TODO: check definition - same node_num/bias, same direction

    # extract blocks boundaries
    boundary_ids = block[:bounds]
    boundaries = bounds[abs.(boundary_ids)]

    # extract 
    points, functions = blockdef2fun(mesh, block, bounds, nodeMapping)

    P12, P14, P34, P32 = points
    c1, c2, c3, c4 = functions
    
    # calculate node positions

    us = getFac(boundaries[1], flipdirection=sign(boundary_ids[1])==-1)
    vs = getFac(boundaries[2], flipdirection=sign(boundary_ids[2])==-1)

    nodearray = calculateTransfiniteNodes(c1, c2, c3, c4, P12, P14, P34, P32, us, vs)

    # add nodes and create elements
    N1 = boundaries[1].node_num
    N2 = boundaries[2].node_num

    last_node_ids = getBoundaryNodeIds(boundary_ids[4], nodeMapping)[end:-1:1]
    node_ids = zeros(Int64, size(last_node_ids, 1))

    for i in range(2, N1-1)
        # add bound node to node ids

        for j in range(2, N2-1)
            nodecoords = collect(nodearray[i, j])
            node_ids[j] = appendNodes!(mesh, [nodecoords,])[1]
        end

        # add bound node to node ids
        node_ids[1] = getBoundaryNodeIds(boundary_ids[1], nodeMapping)[i]
        node_ids[end] = getBoundaryNodeIds(boundary_ids[3], nodeMapping)[i]

        # add elements to mesh
        addElements!(mesh, last_node_ids, node_ids)

        last_node_ids[1:end] = node_ids[1:end]

    end

    node_ids[1:end] = getBoundaryNodeIds(boundary_ids[2], nodeMapping)

    addElements!(mesh, last_node_ids, node_ids)

end
    
function meshTransitionUnit!(mesh::Mesh, nodes1::Vector{Int64}, nodes2::Vector{Int64}, left::Bool)
    cornerpoints_l = mesh.nodes[[nodes1[1], nodes1[3], nodes2[1], nodes2[2]]]
    cornerpoints = mesh.nodes[[nodes1[1], nodes1[end], nodes2[1], nodes2[end]]]
    cornerpoints_r = mesh.nodes[[nodes1[3], nodes1[end], nodes2[2], nodes2[end]]]

    center_ids = zeros(Int64, 3)

    for (i, cp) in enumerate((cornerpoints_l, cornerpoints, cornerpoints_r))
        centerpoint = sum(cp)/size(cp,1)

        center_ids[i] = appendNodes!(mesh, [centerpoint,])[1]
    end

    # add elements

    n1 = nodes1
    n2 = nodes2
    c = center_ids

    appendElement!(mesh, [n1[1], n1[2], c[1], n2[1]])
    appendElement!(mesh, [n1[2], n1[3], c[2], c[1]])
    appendElement!(mesh, [c[1], c[2], n2[2], n2[1]])
    appendElement!(mesh, [n1[3], n1[4], c[3], c[2]])
    appendElement!(mesh, [c[2], c[3], n2[3], n2[2]])
    appendElement!(mesh, [n1[4], n1[5], n2[3], c[3]])
end

function meshTransitionBlock!(mesh::Mesh, block, bounds, nodeMapping)
    
    # extract blocks boundaries
    boundary_ids = block[:bounds]
    boundaries = bounds[abs.(boundary_ids)]

    # extract 
    points, functions = blockdef2fun(mesh, block, bounds, nodeMapping)

    P12, P14, P34, P32 = points
    c1, c2, c3, c4 = functions

    # interpolate border nodes
    N = boundaries[1].node_num
    N2 = boundaries[2].node_num

    # get node positions
    vs = getFac(boundaries[2])

    # storage vector for node_id_rows
    rows = Vector{Vector{Int64}}()

    # left and right node 
    nodes_l = getBoundaryNodeIds(boundary_ids[4], nodeMapping)[end:-1:1]
    nodes_r = getBoundaryNodeIds(boundary_ids[2], nodeMapping)

    println("HallO")
    println([mesh.nodes[n] for n in nodes_l])
    println([mesh.nodes[n] for n in nodes_r])
    println()

    # add first row (boundary)
    push!(rows, getBoundaryNodeIds(boundary_ids[1], nodeMapping))

    # add nodes in rows
    for i in 2:N2-1

        # new number of elements is half of previous
        # number of elements is node_num: N minus 1
        N = Int(floor((N-1)/2))+1

        # interpolation values for transfinite interpolation
        us = collect(1:N-1)[1:end-1]/(N-1)
        vs_cur = [vs[i]]

        # calculate inner node coordinates using transfinite interpolateion
        nodearray = calculateTransfiniteNodes(c1, c2, c3, c4, P12, P14, P34, P32, us, vs_cur)

        # convert nodearray in appropriate from for adding to mesh 
        node_coords = [collect(nodearray[i]) for i in 1:size(nodearray,1)]
        
        # add nodes to mesh
        curnode_ids = appendNodes!(mesh, node_coords)

        # collect node ids of border nodes to get node ids of complete row
        nodeid_l = nodes_l[i]
        nodeid_r = nodes_r[i]

        push!(rows, [nodeid_l, curnode_ids..., nodeid_r])

    end

    push!(rows, getBoundaryNodeIds(boundary_ids[3], nodeMapping)[end:-1:1])


    for row in rows
        println(mesh.nodes[row])
    end

    for i in 2:size(rows,1)
        
        row1 = rows[i-1]
        row2 = rows[i]

        left = true
        println(collect(0:2:size(row2,1)-4))
        for j in 1:(size(row2,1)-1)÷2
            nodes1 = row1[4*(j-1)+1:4*j+1]
            nodes2 = row2[2*(j-1)+1:2*j+1]

            meshTransitionUnit!(mesh, nodes1, nodes2, left)

            left = !left
        end
    end
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
                
