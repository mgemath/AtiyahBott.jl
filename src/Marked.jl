function invert_marks(mark::Marks_type, n_vert::Int64)::Dict{Int64,Vector{Int64}}
    
    D = Dict{Int64,Vector{Int64}}([v for v in 1:n_vert] .=> [Int64[] for _ in 1:n_vert])

    foreach(t -> push!(D[mark[t]], t) , eachindex(mark))

    return D
end

function num_marks(mark::Marks_type, v::Int64)::Int64

    length(mark) == 0 && return 0

    return count(t -> t == v, mark)
end

function ismin(ls::Vector{Int64}, col::Tuple{Vararg{Int64}}, m::Vector{Int64}, par::Vector{Int64}, sub_end::Vector{Int64})::Bool

    marks = sort(unique(m))

    for j in marks # eachindex(marks)
        # find the all the left vertex with the same depth
        for k in findall(i -> (i < j) && ls[i] == ls[j] && col[i] == col[j] && ls[sub_end[i]] == ls[sub_end[j]], eachindex(ls))
            # k = findlast(i -> (i<j) && ls[i] == ls[j], 1:n)

            # if there is no such vertex, then this particular j is minimal
            # k === nothing && continue

            # if both have marks, then this particular j is minimal
            # k in marks && continue

            # necessary conditions for an automorphism
            # ls[sub_end[k]] != ls[sub_end[j]] && continue
            # ls[subgraph_ends[k]] != ls[subgraph_ends[j]] && continue

            # we run backward the graphs to see if there is an isomorphism sending k to j

            # first, we find the last root of k that is not a root of j, and viceversa
            # note that at this point we do not need of k anymore

            par_j::Int64 = copy(j)
            par_k::Int64 = copy(k)
            root_k::Int64 = 0
            root_j::Int64 = 0

            while par_k != par_j
                root_k = par_k
                root_j = par_j

                par_k = par[par_k]
                par_j = par[par_j]
            end

            # Now we can check if there is such isomorphism

            indeces_k::UnitRange{Int64} = root_k:sub_end[root_k]
            indeces_j::UnitRange{Int64} = root_j:sub_end[root_j]

            if view(ls, indeces_k) == view(ls, indeces_j)  # check if the two subtrees are isomorphic
                if col[indeces_k] == col[indeces_j] # check if the colorations are preserved
                    for (x, y) in Base.Iterators.zip(indeces_k, indeces_j)
                        if !(x in marks) && !(y in marks)
                            continue
                        end

                        dim_dif = 0 #num_marks(m,y) - num_marks(m,x)
                        for i in m
                            if i == y
                                dim_dif += 1
                            end
                            if i == x
                                dim_dif -= 1
                            end
                        end

                        if dim_dif > 0
                            return false
                        end

                        # so num_marks(m,y) == num_marks(m,x)
                        # findfirst(i-> i == x, m) is the minimum of all marks pointing at x
                        if dim_dif == 0
                            if findfirst(i -> i == y, m) < findfirst(i -> i == x, m)
                                return false
                            end
                        end
                        break
                    end
                end
            end
        end
        # here means that the isomorphism does not exists
    end
    return true
end