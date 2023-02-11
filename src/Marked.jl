function invert_marks(mark::Tuple{Vararg{Int64}}, n_vert::Int64)::Dict{Int64,Vector{Int64}}
    
    D = Dict{Int64,Vector{Int64}}([v for v in 1:n_vert] .=> [Int64[] for _ in 1:n_vert])

    foreach(t -> push!(D[mark[t]], t) , eachindex(mark))

    return D
end

function num_marks(mark::Tuple{Vararg{Int64}}, v::Int64)::Int64

    length(mark) == 0 && return 0

    return count(t -> t == v, mark)
end