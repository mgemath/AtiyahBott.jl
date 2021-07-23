
const Read = Dict("1"=>1,"2"=>2,"3"=>3,"4"=>4,"5"=>5,"6"=>6,"7"=>7,"8"=>8,"9"=>9,"a"=>10,           
            "b"=>11,"c"=>12,"d"=>13,"e"=>14,"f"=>15,"g"=>16,"h"=>17,"i"=>18,"j"=>19,"k"=>20,
            "l"=>21,"m"=>22,"n"=>23,"o"=>24,"p"=>25,"q"=>26,"r"=>27,"s"=>28,"t"=>29,"u"=>30)
"""
The structure `graph_coloring` define the colorations of a graph.
# Arguments
- `graph::SimpleGraph`: the graph.
- `num_cols::UInt8`: the number of colors.
- `current_color::Array{UInt8,1}`: the array of colors.
"""
mutable struct graph_coloring
    graph::SimpleGraph
    num_cols::UInt8
    current_color::Array{UInt8,1}
end

function graph_coloring( G::SimpleGraph, num_cols::UInt8 )::graph_coloring
    return graph_coloring( G, num_cols, UInt8[] )
end

mutable struct graph_coloring_from_file
    file_name::String
    #color_file::IOStream
    color_file_gz::CodecZlib.TranscodingStreams.TranscodingStream{GzipDecompressor, IOStream}
    current_color::Array{UInt8,1}
    current_aut::Int64
end

# function graph_coloring_from_file( file_name::String )
#     color_file = open( file_name )
#     return graph_coloring_from_file( file_name, color_file, [], 1 )
# end


function graph_coloring_from_file( file_name::String )::graph_coloring_from_file
    color_file_gz = GzipDecompressorStream(open(file_name))
    return graph_coloring_from_file( file_name, color_file_gz, UInt8[], 1 )
end

function Base.iterate( GC::graph_coloring_from_file, c=0 )::Union{Nothing, Tuple{Vector{UInt8}, Int64}}

    st = readline( GC.color_file_gz )
    if st == "STOP" 
        close( GC.color_file_gz )
        return nothing
    end 
    
    st = split( st, "," )
    #GC.current_color = [ parse( UInt8, Read[x] ) for x in split( st[1], "" )]
    GC.current_color = [ UInt8(Read[x]) for x in split( st[1], "" )]
    GC.current_aut = parse( Int64, st[2] )
    return GC.current_color, 0
end

function exists_file_with_colorings(pruf_str::String, n::Int64)::Tuple{Bool,Union{String,Nothing}}
    
    dir = dirname(current_dir)*"/Data/Dimension$n/" #path of the folder containing the colorations

    if isdir( dir ) && pruf_str*"z.gz" in readdir( dir )
        return true, dir*pruf_str*"z.gz"
    else
        return false, nothing
    end
end

function is_coloring(G::SimpleGraph, cols::Array{UInt8,1})::Bool

    num_v = nv( G )
    for i in 1:num_v
        for j in (i+1):num_v
            if cols[i] == cols[j] && has_edge( G, i, j ) 
                return false
            end
        end
    end

    return true
end

function smallest_coloring(G::SimpleGraph, num_cols::UInt8)::Array{UInt8,1}

    cols = [ 0x1 for _ in 1:nv(G)]
    while true
        if is_coloring( G, cols )
            return cols
        end

        next_tuple!( cols, num_cols )
    end
end


function Base.iterate(GC::graph_coloring, c=0)::Union{Nothing, Tuple{Vector{UInt8}, Int64}}

    if GC.current_color == []
        GC.current_color = smallest_coloring( GC.graph, GC.num_cols )
        return GC.current_color, 0
    end

    while true
        v = next_tuple!( GC.current_color, GC.num_cols )
        if v == false
            return nothing
        end

        if is_coloring( GC.graph, GC.current_color )
            return GC.current_color, 0
        end

    end
end
    

function next_tuple!( tuple::Array{UInt8,1}, max_entry::UInt8 )::Union{Nothing, Bool}

    k = length( tuple )
    for k in k:-1:0

        if k == 0  
            return false        
        elseif tuple[k] < max_entry 
            tuple[k] += 1
            break
        else
            tuple[k] = 1
        end
    end

    return nothing
end


"""
Compute the graph with a specific Prufer sequence.
"""
function PruferToGraph(prufer::Vector{UInt8})::SimpleGraph
    
    num_v = length(prufer)+2 #number of vertices
    g = Graph(num_v) #initialize a graph with specific number of vertices
    
    degree = [count(i->(i==j),prufer) for j in 1:num_v ] #new version of degree

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

"""
Extract a Prufer sequence and a number from a string. The number is the number of automorphisms of the graph with that Prufer sequence.
"""
function get_graph(str::String)::Tuple{SimpleGraph{Int64}, Int64}

    s = split(str, ',')
    g = PruferToGraph(UInt8[ Read[string(s[1][i])] for i in 1:length(s[1])]) #get the graph
    a = parse(Int64,s[2])       
    return (g, a)
end

"""
Return all the arrays of length `l` such that the sum of all elements of the array is `d`.
"""
function get_weights(l::Int64, d::Int64)::Vector{Vector{Int64}}
    return vcat(unique.(permutations.(partitions(d,l)))...)
end