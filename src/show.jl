using Plots

function show!(p::Plots.Plot, boundary::StraightBoundary, meshdef::MeshDef)
    vert_start = boundary.vertice_start
    vert_end = boundary.vertice_end
    coords = meshdef.vertices[[vert_start, vert_end]]

    plot!(p, [coords[1][1], coords[2][1]], [coords[1][2], coords[2][2]], linecolor="blue")
end

function show!(p::Plots.Plot, boundary::CircularBoundary, meshdef::MeshDef) 
    vert_start = boundary.vertice_start
    vert_end = boundary.vertice_end
    coords = meshdef.vertices[[vert_start, vert_end]]
    
    # two circle intersection -> https://math.stackexchange.com/a/1367732
    dist = norm(coords[1]-coords[2])
    midpoint = (coords[1]+coords[2])/2
    normaldist = sqrt(boundary.radius^2 - dist^2/4)
    normalvec = reverse((coords[1] - coords[2])) .* [1 -1]
    center = midpoint .+ sign(boundary.radius) .* normalvec .* normaldist 
    
    leg1 = coords[1] .- center
    leg2 = coords[2] .- center
    
    γ = angle(leg1, leg2)
    
    x(α) = leg1[1] * cos(α) - leg1[2] * sin(α)
    y(α) = leg1[1] * sin(α) + leg1[2] * cos(α)
    
    plot!(p, x, y, 0, γ, leg=false)
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