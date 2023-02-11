struct ColorsIt
    ls::Vector{Int64}
    max_col::Int64
end

function Base.iterate( CI::ColorsIt, state::Tuple = () )::Union{Nothing, Tuple{NTuple{length(CI.ls), UInt8}, Tuple}}

    if isempty(state)
        c::Vector{UInt8} = my_minimal_coloring(CI.ls)
        # return ntuple(i -> c[i], length(c)), (init_ColorsIt(CI.ls), c)
        return (c..., ), (init_ColorsIt(CI.ls), c)
    end

    ((rev_dfs, parents, has_ci, subgraph_ends_rev, subgraph_ends, left_siblings, n), col) = state

    ls = CI.ls
    max_col = CI.max_col

    # col::Vector{Int64} = copy(state)

    # has_ci::Bool = get_prop( g, :has_ci )
    
    # # set up the maximum color, number of vertices, and the parents for the vertices
    # n::Int64 = length( ls )
    # parents::Vector{Int64} = get_prop( g, :parents )

    # # set up the reverse dfs order for the vertices, the end of subgtrees originating from each vertex and 
    # # the left siblings of each vertices
    # rev_dfs::Vector{Int64} = get_prop( g, :reverse_dfs )
    # subgraph_ends::Vector{Int64} = get_prop( g, :ends_of_subgraphs )
    # subgraph_ends_rev::Vector{Int64} = get_prop( g, :ends_of_subgraphs_rev )
    # left_siblings::Vector{Int64} = get_prop( g, :left_siblings )

    
    #=
        A vertex is colorable if it hasn't reached its maximum possible color.
            
        The vertex 1 is colorable if either 
        1.) The graph has no central involution and the color is 1 is not the max color; or 
        2.) The graph has central involutiuon and its color is not the second largest color. 
    =#
    
    # this variable shows how far we need to look in the tree to find the vertex whose 
    # color is modified
    look_end::Int64 = subgraph_ends[1]

    # in colorable we keep the index of the vertex that is to be colored
    colorable::Int64 = 0 #colorable_1 ? 1 : 0

    # we go through the vertices in reverse DFS order
    k::Int64 = 2

    while k <= n
        v::Int64 = rev_dfs[k]
        #println( k )

        # compute end point of the subgraph T_v stemming at vertex v
        v_end::Int64 = subgraph_ends[v]
        v_end_rev::Int64 = subgraph_ends_rev[v]
        # compute the root left sibling of the subgraph T_v
        # (the subgraph stemming from a vertex with the same parent as v, left of v)
        left_s::Int64 = left_siblings[v]
        
        left_sibling::UnitRange{Int64} = 1:0

        if left_s != -1
            # the left sibling exists
            # we compute the range of indices that correspond to the left sibling
            left_sibling = left_s:(v-1)
        end 

        # check if T_v is isomorphic to the left sibling and if it has the same coloring
        # If yes, then the coloring of T_v cannot be increased and hence this part of 
        # the tree can be ignored

        if view(ls,v:v_end) == view(ls,left_sibling) && view(col,v:v_end) == view(col,left_sibling)
            #println( k )
            k = findfirst( x -> rev_dfs[x] == left_s, 1:n )
            #println( v, left_s )
            #println( "jumping to $k" )
            continue 
        end 

        # we decide if the color of v can be increased
        #println( "v is $v")
        if col[v] < max_col - 1 
            # we update colorable if necessary 
            if colorable < v colorable = v end 
            #we update look_end 
            look_end = v_end_rev 
        elseif col[v] == max_col - 1 && col[parents[v]] != max_col
            if colorable < v colorable = v end
            look_end = v_end_rev 
        end 
        
        if v == look_end && colorable != 0
            # we got the end of a branch that has colorable vertex
            # we need to look no further
            break 
        end  
        k += 1
    end 

    # check if one is colorable 
    colorable_1::Bool = has_ci ? col[1] < max_col - 1 : col[1] < max_col

    # if did not find colorable vertex then return nothing
    if colorable == 0 
        if colorable_1 
            colorable = 1
        else 
            return nothing
        end
    end

    # if colorable == 0 && !colorable_1 
    #     return nothing 
    # elseif colorable == 0 
    #     colorable = 1
    # end
    
    if colorable != 1 && col[parents[colorable]] == col[colorable] + 1
        # the parent of v already has color col[v] + 1 and so we increase color by two
        col[colorable] += 0x02 
    else 
        # else we increase the color by one
        col[colorable] += 0x01
    end 

    # we reset the coloring of each subtree on the right side of v
    # starting from v+1

    j::Int64 = colorable+1
    while j <= n
        root_color::UInt8 = 0x00

        # find the end of the subtreee T_j stemming from vertex j
        es::Int64 = subgraph_ends[j]

        # determine the color for the root j of this subtree
        if has_ci && j == 2 
            # if has central involution and j is vertex two, then its color
            # is set to the color of vertex 1 plus 1 
            root_color = col[1] + 0x01
        # else 
        #     # else no constaint on root color
        #     root_color = 0x00
        end 

        # compute the minimal coloring for the subtree T_j
        # with parent color being the color of parent[j]
        col[j:es] = my_minimal_coloring( ls[j:es], parent_color = col[parents[j]], root_color = root_color )

        # the following subtree that needs to be dealt with stems from the end of T_j plus 1 
        j = es+1
    end

    new_col::NTuple{length(ls), UInt8} = (col..., ) # ::NTuple{length(ls), UInt8} = ntuple(i -> col[i], length(ls))

    # return the coloring computed
    return new_col, ((rev_dfs, parents, has_ci, subgraph_ends_rev, subgraph_ends, left_siblings, n), col)# (g, col)
end 

#########AUX FUNCTIONS#######

function init_ColorsIt(ls::Vector{Int64})::Tuple

    n::Int64 = length(ls)
    
    par::Vector{Int64} = [0 for _ in 1:n]

    foreach(v -> par[v] = findlast(i -> i < v && ls[i] == ls[v] - 1, eachindex(ls)), 2:n)

    stack::Vector{Int64} = [1]
    rev_dfs::Vector{Int64} = []
    v::Int64 = 0
    while length(stack) > 0 
        v = pop!( stack )
        push!( rev_dfs, v )
        append!( stack, findall(x-> par[x] == v, 1:n))
    end 

    ans = (
            rev_dfs, # reverse_dfs
            par, # parents
            my_has_central_involution( ls ), # has_ci
            [ end_of_subgraph_rev( ls, rev_dfs, x ) for x in 1:n ], # ends_of_subgraphs_rev
            [ my_end_of_subgraph( ls, x ) for x in 1:n ], # ends_of_subgraphs
            [ my_root_of_left_sibling_subtree( par, x ) for x in 1:n ], # left_siblings    
            n # length( seq )
    )

    return ans
end

function my_minimal_coloring( ls::Vector{Int64}; parent_color::UInt8 = 0x00, root_color::UInt8 = 0x00 )::Vector{UInt8}#::NTuple{length(ls), UInt8}
    
    # choose the root color
    # col_pair::NTuple{2, UInt8}
    
    if root_color == 0x00
        if parent_color == 0x00
            root_color = 0x01
        else
            root_color = parent_color == 0x01 ? 0x02 : 0x01 
        end
    end

    if (root_color == 0x01) == isodd(ls[1])
        col_pair = (0x02, 0x01) # = col_even, col_odd
    else 
        col_pair = (0x01, 0x02)
    end

    ans::Vector{UInt8} = [col_pair[(ls[i] % 2) + 1] for i in eachindex(ls)]
    # for NTuple use ans::NTuple{length(ls),UInt8} = col_pair[(ls .% 2) .+ 1] 
    
    ans[1] = root_color

    return ans
end 

function my_root_of_left_sibling_subtree( par::Vector{Int64}, v::Int64  )::Int64
    
    p::Int64 = par[v]
    if p == 0
        return -1
    end 
    desc::Vector{Int64} = [ x for x in my_children_vertices( par, p ) if x < v ]
    
    if isempty(desc)
        return -1
    end 

    return maximum( desc )
end 

function my_children_vertices( par::Vector{Int64}, v::Int64 )::Vector{Int64}

    return findall(i -> par[i] == v, eachindex(par))
end

function my_end_of_subgraph( ls::Vector{Int64}, v::Int64 )

    # find the first vertex whose level is not grater then the level of v
    # if found nothing, then return the length of ls
    # otherwise return end_ver+v-1
    return findfirst( i -> i == length(ls) || (i >= v && ls[i+1] <= ls[v]), eachindex(ls))
end 

function my_end_of_subgraph_rev( ls::Vector{Int64}, r_dfs::Vector{Int64}, v::Int64 )::Int64

    # position of v in r_dfs 
    ps = findfirst( x -> r_dfs[x] == v, 1:length( ls ))
    end_ver::Union{Nothing,Int64} = findfirst( k -> ls[k] <= ls[v], view(r_dfs,ps+1:length(r_dfs)) )

    # if found nothing, then return the length of ls
    # otherwise return end_ver+v-1
    return end_ver === nothing ? r_dfs[end] : r_dfs[ps+end_ver-1]
end 

function end_of_subgraph_rev( ls::Vector{Int64}, r_dfs::Vector{Int64}, v::Int64 )::Int64 #to review

    # position of v in r_dfs 
    ps = findfirst( x -> r_dfs[x] == v, 1:length( ls ))
    end_ver::Union{Nothing,Int64} = findfirst( k -> ls[k] <= ls[v], r_dfs[ps+1:end] )

    # if found nothing, then return the length of ls
    # otherwise return end_ver+v-1
    return end_ver === nothing ? r_dfs[end] : r_dfs[ps+end_ver-1]
end 

# function my_end_of_subgraph_rev( ls::Vector{Int64}, r_dfs::Vector{Int64}, v::Int64 )::Int64

#     # position of v in r_dfs 
#     ps = findfirst( x -> r_dfs[x] == v, 1:length( ls ))
#     k = findfirst( i -> i == length(ls) || (i >= ps && ls[i+1] <= ls[v]), eachindex(r_dfs))
#     return r_dfs[k]
# end 

function my_has_central_involution( ls::Vector{Int64} )::Bool

    # if the graph is o--o then the answer is yes
    if ls == [1,2] return true end
    
    if iseven(length( ls ))
        two = findfirst(i -> i > 2 && ls[i] == 2, eachindex(ls))
        if view(ls,3:two-1) == 1 .+ view(ls,two:length(ls))
            return true
        end
    end

    return false
end

######Of Iterators#########

function Base.eltype(CI::ColorsIt)
    NTuple{length(CI.ls), UInt8}
end

function Base.length(CI::ColorsIt)::Int64
    function n_colorations(ls::Vector{Int64}, n::Int64)::Int64

        if length(ls) < 3
            return ls[1] == 1 ? binomial(n, length(ls)) : (n-1)^length(ls) # ls[1] == 1 is equivalent to be a starting-point graph
        end
    
        my_child = findall(i -> i>length(ls) || ls[i]==ls[1]+ 1, 1:length(ls)+1)  # find all the children of the root, plus length(ls)+1
        
    
        if ls[1] == 1  # ls[1] == 1 is equivalent to be a starting-point graph
            if iseven(length(ls))  # necessary condition to have the bad involution
                if view(ls,my_child[1]+1:my_child[2]-1) == 1 .+ view(ls,my_child[2]:length(ls))  # check if it has the bad involution
                    c = n_colorations(ls[my_child[1]:my_child[2]-1], n)   # compute the number of colorations of the main subgraph
                    return div(n*(c^2),2*(n-1))  # return the number of coloration computing only the number of colorations of the main subtree
                end
            end
            ans = n  # in the starting-point graph, the root can assume n values since it has no parent
        else
            ans = n - 1 # we are not in the starting-point graph, so the root has a parent to deal with
        end
    
        # here we compute the number of colorations of each subgraph
        last_sub = ls[my_child[1]:my_child[2]-1]  # this is the subgraph
        m = 1                                     # this is its multiplicity
    
        for x in 3:length(my_child)
            if last_sub == view(ls, my_child[x-1]:my_child[x]-1)  # we run up to find a different subgraph
                m += 1
            else
                c = n_colorations(last_sub, n)   # number of colorations of this subgraph...
                ans *= binomial(c+m-1,m)         # ...counted with multiplicity
                last_sub = ls[my_child[x-1]:my_child[x]-1]  # pass to the next subgraph
                m = 1
            end
        end
    
        c = n_colorations(last_sub, n)  # compute the number of colorations of the last subgraph
        ans *= binomial(c+m-1,m)        # with multiplicity
    
        return ans
    end
    return n_colorations(CI.ls, CI.max_col)
end

### this function counts the isomorphisms of a tree with level sequence ls. Optionally, the tree can be colored with coloration col

# function count_iso(ls::Vector{Int64}, col::Vector{UInt8} = UInt8[])::Int64
function count_iso(ls::Vector{Int64}, col::Tuple{Vararg{UInt8}} = ())::Int64

    is_empty::Bool = isempty(col)

    if length(ls) < 3
        if is_empty
            return ls[1] == 1 ? length(ls) : 1 # ls[1] == 1 is equivalent to be a starting-point graph
        else
            return 1
        end
    end

    my_child = findall(i -> i>length(ls) || ls[i]==ls[1]+ 1, 1:length(ls)+1)  # find all the children of the root, plus length(ls)+1
    

    if is_empty    
        if ls[1] == 1 # ls[1] == 1 is equivalent to be a starting-point graph
            if iseven(length(ls))  # necessary condition to have the bad involution
                if view(ls,my_child[1]+1:my_child[2]-1) == 1 .+ view(ls,my_child[2]:length(ls))  # check if it has the bad involution
                    a = count_iso(ls[my_child[1]:my_child[2]-1])   # compute the number of colorations of the main subgraph
                    return 2*(a^2)
                end
            end
        end
    end

    # here we compute the number of colorations of each subgraph
    last_sub::Vector{Int64} = ls[my_child[1]:my_child[2]-1]  # this is the subgraph
    last_col::Tuple{Vararg{UInt8}} = col
	
	if !is_empty
        last_col = col[my_child[1]:my_child[2]-1]
    end
    # last_col::Vector{Int64} = col[my_child[1]:my_child[2]-1]
    m::Int64 = 1
    ans::Int64 = 1
    

    for x in 3:length(my_child)
        if last_sub == view(ls, my_child[x-1]:my_child[x]-1) && (is_empty || last_col == col[my_child[x-1]:my_child[x]-1]) #view(col, my_child[x-1]:my_child[x]-1))
            m += 1
        else
            a = count_iso(last_sub, last_col)   # number of colorations of this subgraph...
            ans *= (a^m)*factorial(m)
            last_sub = ls[my_child[x-1]:my_child[x]-1]  # pass to the next subgraph
            # last_col = col[my_child[x-1]:my_child[x]-1]  # pass to the next subgraph
            if !is_empty
                last_col = col[my_child[x-1]:my_child[x]-1]
            end
            m = 1
        end
    end

    a = count_iso(last_sub, last_col)  # compute the number of colorations of the last subgraph
    ans *= (a^m)*factorial(m)       # with multiplicity

    return ans
end