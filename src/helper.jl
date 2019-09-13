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