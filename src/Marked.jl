mutable struct Marks
    n_vert::Int64
    m::Int64
    get_vertex::Dict{Int64,Int64}
end

function Marks(v::Int64,m::Int64)::Marks

    return Marks(v, m , Dict{Int64, Int64}())
end

function Base.iterate(mark::Marks, c::Int64 = 0)::Union{Nothing, Tuple{Marks, Int64}}
    
    if mark.m == 0 && c == 1 
        return nothing
    end

    if isempty(mark.get_vertex)
        mark.get_vertex = Dict([i for i in 1:mark.m].=> 1)
        return mark, 1
    end

    index = findfirst(x -> mark.get_vertex[x] < mark.n_vert, 1:mark.m)
    if index === nothing
        return nothing
    end

    next_mark!(mark, index)
    return mark, 1
end


function next_mark!(mark::Marks, index::Int64)::Nothing
    
    mark.get_vertex[index] += 1
    for i in 1:(index-1)
        mark.get_vertex[i] = 1
    end
    
    return nothing
end            
#=
"""
Take the marks of a graph and return a dictionary Dict{Int64,Vector{Int64}} which associate to each vertex the set of marks on the vertex.
"""=#
function invert_marks(mark::Marks)::Dict{Int64,Vector{Int64}}
    new_D = Dict{Int64,Vector{Int64}}()
    D = mark.get_vertex
    for k in keys(D)
        D[k] in keys(new_D) ? push!(new_D[D[k]],k) : new_D[D[k]] = [k]
    end
    for v in 1:mark.n_vert
        if !(v in keys(new_D))
            new_D[v]=Int64[]
        end
    end
    return new_D
end
#=
"""
Return the number of marks assigned to some vertex v.
# Arguments
- `mark::Marks`: the marks assigned to a graph.
- `v::Int64`: the vertex.
"""=#
function num_marks(mark::Marks, v::Int64)::Int64
    
    mark.m == 0 && return 0

    return count(t -> mark.get_vertex[t] == v, 1:mark.m)
end