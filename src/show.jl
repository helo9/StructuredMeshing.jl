using Plots

function show!(p::Plots.Plot, boundary::Boundary, meshdef::MeshDef)
    vert_start = boundary.vertice_start
    vert_end = boundary.vertice_end
    coords = meshdef.vertices[[vert_start, vert_end]]

    plot!(p, [coords[1][1], coords[2][1]], [coords[1][2], coords[2][2]], linecolor="blue")
end


function show!(p::Plots.Plot, vertice::Vector{Float64})
    scatter!(p, [vertice[1]], [vertice[2]], marker=(:circle, :red, 6))
end

function show!(p::Plots.Plot, element::Vector{Int64}, mesh::Mesh)
    start_node = mesh.nodes[element[1]]
    
    for node_id in element[2:end]
        
        cur_node = mesh.nodes[node_id]
        
        plot!(p, [start_node[1], cur_node[1]], [start_node[2], cur_node[2]], linecolor="grey")
        
        start_node = cur_node
    end
    
    last_node = mesh.nodes[element[1]]
    plot!(p, [start_node[1], last_node[1]], [start_node[2], last_node[2]], linecolor="grey")
    
end

function show!(p::Plots.Plot, mesh::Mesh)
    
    # TODO vectorize!
    for node in mesh.nodes
        scatter!(p, [node[1]], [node[2]], marker=(:circle, :blue, 3))
    end
    
    for element in mesh.elements
        show!(p, element, mesh)
    end
end

function show(meshdef::MeshDef; hideVertices::Bool=false)
    p = plot(legend=false)
    
    for boundary in meshdef.bounds
        show!(p, boundary, meshdef)
    end
    
    mesh = emptyMesh()
    meshBoundaries!(mesh, meshdef)
    show!(p, mesh)
    
    if !hideVertices
        for vertice in meshdef.vertices
            show!(p, vertice)
        end
    end
    
    p
end

function show(mesh::Mesh)
    p = plot(legend=false)
    
    show!(p, mesh)
    
    p
end