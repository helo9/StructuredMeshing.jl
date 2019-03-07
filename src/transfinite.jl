function calculateTransfiniteNodes(c1::Function, c2::Function, c3::Function, c4::Function, P12::Vector{Float64}, 
    P14::Vector{Float64}, P34::Vector{Float64}, P32::Vector{Float64}, us::Vector{Float64}, vs::Vector{Float64})

S = Array{Union{Nothing, Tuple}}(nothing, size(us,1), size(vs,1))

for (i, u) in enumerate(us)
    for (j, v) in enumerate(vs)
        tmp = (1-v)*c1(u) + v*c3(u) + (1-u)*c2(v) + u*c4(v) - ( (1-u)*(1-v)*P12 + u*v*P34 + u*(1-v)*P14 + (1-u)*v*P32 )
        S[i,j] = tuple(tmp...)
    end
end

return S

end