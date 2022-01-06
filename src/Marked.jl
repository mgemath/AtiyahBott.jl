struct Marks
    n_vert::Int64
    m::Int64
end
Base.eltype(::Type{Marks}) = Vector{Int}
Base.length(mark::Marks) = mark.n_vert^mark.m

#each element of the cycle for is an array of length m such that the position at i is the vector containing the mark I

function Base.iterate(mark::Marks, prev::Union{Nothing, Vector{Int}} = nothing)::Union{Nothing, Tuple{Vector{Int}, Vector{Int}}}
    
    prev === nothing && return (ones(Int,mark.m), ones(Int,mark.m))

    index = findfirst(x -> prev[x] < mark.n_vert, 1:mark.m)
    
    index === nothing && return nothing

    prev[index] += 1
    for i in 1:(index-1)
        @inbounds prev[i] = 1
    end

    return (prev, prev)
end

function invert_marks(mark::Vector{Int}, n_vert::Int)::Dict{Int64,Vector{Int64}}
    
    D = Dict{Int64,Vector{Int64}}([v for v in 1:n_vert] .=> [Int[] for _ in 1:n_vert])

    for t in 1:length(mark)
        push!(D[mark[t]], t)
    end
    return D
end

function num_marks(mark::Vector{Int}, v::Int64)::Int64

    length(mark) == 0 && return 0

    return count(t -> mark[t] == v, 1:length(mark))
end