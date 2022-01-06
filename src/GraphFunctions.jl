struct graph_coloring
    graph::SimpleGraph
    num_cols::UInt8
end

Base.eltype(::Type{graph_coloring}) = Vector{UInt8}
Base.length(GC::graph_coloring) = GC.num_cols*((GC.num_cols-1)^(nv(GC.graph)-1))

function Base.iterate(GC::graph_coloring, prev::Vector{UInt8} = UInt8[])::Union{Nothing, Tuple{Vector{UInt8}, Vector{UInt8}}}

    if isempty(prev)
        prev = ones(UInt8, nv(GC.graph))
    end

    while true

        index = findfirst(x -> prev[x] < GC.num_cols, 1:length(prev))
    
        index === nothing && return nothing

        prev[index] += UInt8(1)
        for i in 1:(index-1)
            @inbounds prev[i] = UInt8(1)
        end
        
        #is_coloring(GC.graph, prev) && break
        all(e -> prev[src(e)] != prev[dst(e)], edges(GC.graph)) && break
    end

    return prev, prev #use copy(prev), prev in order to collect the elements
end

########

const Read = Dict("1"=>1,"2"=>2,"3"=>3,"4"=>4,"5"=>5,"6"=>6,"7"=>7,"8"=>8,"9"=>9,"a"=>10,           
            "b"=>11,"c"=>12,"d"=>13,"e"=>14,"f"=>15,"g"=>16,"h"=>17,"i"=>18,"j"=>19,"k"=>20,
            "l"=>21,"m"=>22,"n"=>23,"o"=>24,"p"=>25,"q"=>26,"r"=>27,"s"=>28,"t"=>29,"u"=>30)

struct graph_coloring_from_file
    color_file_gz::CodecZlib.TranscodingStreams.TranscodingStream{GzipDecompressor, IOStream}
end

Base.eltype(::Type{graph_coloring_from_file}) = Tuple{Vector{UInt8}, Int64}
# function Base.length(GC::graph_coloring_from_file)
#     res = countlines(GC.color_file_gz) - 1
#     seekstart(GC.color_file_gz)
#     return res
# end
function Base.iterate(GC::graph_coloring_from_file, ::Int64 = 0)::Union{Nothing, Tuple{Tuple{Vector{UInt8}, Int64}, Int64}}

    st = readline(GC.color_file_gz)
    if st == "STOP" 
        close(GC.color_file_gz)
        return nothing
    end 
    
    st_s = split(st, ",")
    
    current_color = [UInt8(Read[x]) for x in split(st_s[1], "")]
    current_aut = parse(Int64, st_s[2])

    return (current_color, current_aut), 0
end

#=
"""
Compute the graph with a specific Prufer sequence.
"""=#
function PruferToGraph(prufer::Vector{UInt8})::SimpleGraph
    
    num_v = length(prufer)+2 #number of vertices
    g = Graph(num_v) #initialize a graph with specific number of vertices
    
    degree::Vector{Int64} = zeros(Int64,num_v)

    for i in 1:num_v-2
        degree[prufer[i]] += 1
    end

    for i in 1:(num_v-2)
        for j in 1:num_v
            if degree[j] == 0 #If j is not present in prufer set 
                add_edge!(g, j, prufer[i]) #add edge
                degree[j] = -1 #Remove from Prufer
                degree[prufer[i]] -= 1 #low the degree of the vertex i
                break
            end
        end
    end

    # For the last element, we look at those vertices still with degree 
    i = indexin(0, degree)[1]
    j = num_v + 1 - indexin(0, reverse(degree))[1]
    add_edge!(g, i, j)

    return g
end
#=
"""
Extract a Prufer sequence and a number from a string. The number is the number of automorphisms of the graph with that Prufer sequence.
"""=#
# function get_graph(str::String)::Tuple{SimpleGraph{Int64}, Int64}

#     s = split(str, ',')
#     g = PruferToGraph(UInt8[Read[string(s[1][i])] for i in 1:length(s[1])]) #get the graph
#     a = parse(Int64,s[2]) 

#     return (g, a)
# end
#=
"""
Return all the arrays of length `l` such that the sum of all elements of the array is `d`.
"""=#
function get_weights(l::Int64, d::Int64)::Vector{Vector{Int64}}

    return vcat(unique.(permutations.(partitions(d,l)))...)
end

struct graph_gen
    list_g::IOStream
    tot_number::Int64
    in_dir::Union{Vector{String},Nothing}
end

Base.eltype(::Type{graph_gen}) = Tuple{String, Int64, Bool}
Base.length(Gg::graph_gen) = Gg.tot_number

function Base.iterate(Gg::graph_gen, this_num::Int64 = 1)::Union{Nothing, Tuple{Tuple{String, Int64, Bool}, Int64}}
    
    if Gg.tot_number < this_num 
        close(Gg.list_g)
        return nothing
    end

    this_num += 1

    st = split(readline(Gg.list_g), ",")

    has_file::Bool = true

    if Gg.in_dir !== nothing
        has_file = st[1]*"0.gz" in Gg.in_dir
    end

    top_aut::Int64 = parse(Int64, st[2])

    return (string(st[1]), top_aut, has_file), this_num
end