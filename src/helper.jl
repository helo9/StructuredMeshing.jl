function circlecenter(coords1, coords2, radius)
    # two circle intersection -> https://math.stackexchange.com/a/1367732
    dist = norm(coords1-coords2)
    midpoint = (coords1+coords2)/2
    normaldist = sqrt(radius^2 - dist^2/4)
    normalvec = reverse((coords1 - coords2)) .* [1, -1]
    normalvec /= norm(normalvec)
    center = midpoint .+ sign(radius) .* normalvec .* normaldist
    
    return center
end

function circlesection(coords1, coords2, center, u)
    leg1 = coords1 .- center
    leg2 = coords2 .- center
    
    γ = angle(vec(leg1), vec(leg2))
    
    x(au) = leg1[1] * cos(au*γ) - leg1[2] * sin(au*γ) + center[1]
    y(au) = leg1[1] * sin(au*γ) + leg1[2] * cos(au*γ) + center[2]
    
    return (x.(u), y.(u))
end

function boundaryfun(boundary::StraightBoundary, vertices, inverse_direction=false)
    coords1 = vertices[boundary.vertice_start]
    coords2 = vertices[boundary.vertice_end]
    
    connection = coords2 - coords1
    
    if inverse_direction
        return coords2, u -> coords2 - connection .* u
    else
        return coords1, u -> coords1 + connection .* u
    end
end

function boundaryfun(boundary::CircularBoundary, vertices, inverse_direction=false)
    coords1 = vertices[boundary.vertice_start]
    coords2 = vertices[boundary.vertice_end]
    
    center = circlecenter(coords1, coords2, boundary.radius)
    
    leg1 = coords1 .- center
    leg2 = coords2 .- center
    
    γ = angle(vec(leg1), vec(leg2))
    
    T(u) = [cos(u*γ) -sin(u*γ); sin(u*γ) cos(u*γ)]
    
    if inverse_direction
        return coords2, u -> T(1-u) * leg1 + center
    else
        return coords1, u -> T(u) * leg1 + center
    end
end