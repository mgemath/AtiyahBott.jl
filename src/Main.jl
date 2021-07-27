export
    AtiyahBottFormula,
    AtiyahBottFormulaForGraph

const  number_trees = [1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]
#the number of non-isomorphic graphs with given number of vertices (starting from 2)

# ProgressData contains the data necessary to keep the progress bar up-to-date
mutable struct ProgressData
    progress_bar::Progress
    current_graph::Int64
    top_aut::Int64
    threshold::Int64
end

# create trivial progress data in case the user does not want progress bar
function EmptyProgressData()::ProgressData
    P::Progress = Progress(0)
    P.enabled = false
    return ProgressData(P, 0, 0, 0)
end

# the following function performs the computation of the Atiyah-Bott formula
# for a particular graph.
"""
    AtiyahBottFormulaForGraph(g, pruf_str, aut, n, deg, n_marks, P, s)

Apply the Atiyah-Bott formula to all colorations of a specific graph. It is useful for splitting the computation in multiple parts, to be computed in single threads.
"""
function AtiyahBottFormulaForGraph(g::SimpleGraph, pruf_str::String, 
    aut::Int64, n::Int64, deg::Int64, n_marks::Int64, P, s::Vector{Rational{BigInt}},
    progress_data::ProgressData = EmptyProgressData())::Vector{Rational{BigInt}}

    local weights::Vector{Vector{Int64}} = get_weights(nv(g)-1, deg) #the array of array of weights

    local n_results::Int64 = length(P)   #this store the number of final results of our computation
    local result::Vector{Rational{BigInt}} = [Rational{BigInt}(0) for _ in 1:n_results] #the array of final results
    
    #for N in 1:n

        from_file, file_name = exists_file_with_colorings( pruf_str, n ) #n <- N

        if from_file 
            cols = graph_coloring_from_file( file_name )        
        else 
            #println( "graph not found!!!" )
            cols = graph_coloring( g, UInt8(n+1) )   #n <- N
        end 


        for c in cols   #we run among all colorations of g
            if from_file
                aut = cols.current_aut  #we are read the number of automorphisms of the colorated graph from the file
            #else
                #N+1 in cols.current_color || continue
            end
            
            for m in marks(nv(g), n_marks)    #we run among all marks of g, if n_marks==0 we have only the empty mark
                for w in weights          #we run among all weights of g
                    try
                        local Euler::Rational{BigInt} = Euler_inv(g,c,w,s,m)//(aut*prod(w)) #the contribuition of the Euler class in the Atiyah-Bott formula
                        for res in 1:n_results      #compute each term of the array P
                            result[res] += P[res](g,c,w,s,m)*Euler    #apply Atiyah-Bott
                        end
                    catch err 
                        if isa(err, DivideError) 
                            error("Some division by zero occurred. Try again")
                        end
                        println(err)
                        error("Some error occurred")
                        return [Rational{BigInt}(0) for _ in 1:n_results]
                    end
                end

                if progress_data.progress_bar.enabled
                
                    progress_data.current_graph += progress_data.top_aut÷aut   
                    #update the progress bar
                    update!(progress_data.progress_bar, 
                            progress_data.current_graph,
                            showvalues = [(:"Total number of graphs",progress_data.threshold),
                            (:"Current graph",progress_data.current_graph)])
                end
            end
        end
    #end
    return result
end 


"""
    AtiyahBottFormula(n, d, m, P, do_check, show_bar, down_col)

Apply the Atiyah-Bott residue formula to the class `P`, in the moduli space of rational marked stable maps to the projective space of dimension `n` of degree `d` with `m` marks.
# Arguments
- `n::Int64`: the dimension of the projective space, it must be between 1 and 254.
- `d::Int64`: the degree of the stable maps, it must be between 1 and 13.
- `m::Int64`: the number of marks.
- `P`: the equivariant class.
- `do_check::Bool`: if `true`, checks if `P` is a well defined zero cycle, and stops the computation if this is not true. If `false`, the computation may have an unexpected behaviour. By default is `true`.
- `show_bar::Bool`: hide the progress bar if and only if this condition is `false`. By default is `true`.
- `down_col::Bool`: check if all colorations needed in the computation are in the folder Data, and download them otherwise. Once downloaded, those files can be used for future computations. By default is `true`.

The general construction of `P` is the following:
```julia-repl
julia> P = (g,c,w,s,m) ->
```
After `->`, one has to write an expression in the equivariant classes. All such equivariant classes are functions starting with `(g,c,w,s)` or `(g,c,w,s,m)`. At the end, they can have more arguments. The expression is a polynomial combination of the equivariant classes. We compute the degree of `P` by

```julia-repl
julia> AtiyahBottFormula(n,d,m,P);
```

# Example
```julia-repl
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);
julia> AtiyahBottFormula(3,1,0,P);
Warning: the class is not a 0-cycle.
julia> AtiyahBottFormula(4,1,0,P);
Result: 2875//1
julia> AtiyahBottFormula(4,1,0,P,false);             #same as before, but without the preliminary check on `P`
julia> AtiyahBottFormula(4,1,0,P,false,false);       #same as before, but without showing the progress bar
julia> AtiyahBottFormula(4,1,0,P,false,false,false); #same as before, but without checking for the colorations
```

The function returns an array of the same dimension of `P` (non-vectorized classes are assumed as 1-dimensional arrays). The Julia notation for accessing to array is `name_of_array[i]` where `i` is an index starting from 1.

# Example
```julia-repl
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)*Hypersurface(g,c,w,s,3);
julia> x = AtiyahBottFormula(3,2,0,P)[1];
Result: 81//1
julia> x
81//1
```

The class `P` supports parameters.
```julia-repl
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,3)*(Incidency(g,c,w,s,2)//3)^(d-1);
julia> d = 2;
julia> AtiyahBottFormula(3,d,0,P);
Result: 27//1
julia> d = 3;
julia> AtiyahBottFormula(3,d,0,P);
Result: 84//1
```

More examples are available in the support of the equivariant classes. It is enough to type `?` and then the name of the class. Currently, the supported classes are:

* `O1_i`         (Euler class of ``ev^*O(1)``, where ``ev`` is a specific evaluation map)
* `O1`           (product of all `O1_i`)
* `Psi`          (cycle of ``psi``-classes)
* `Jet`          (Euler class of the jet bundle ``J^p`` of ``ev^*O(q)``)
* `Hypersurface` (Euler class of the direct image of ``ev^*O(b)``)
* `Incidency`    (cycle parameterizing curves meeting a linear subspace)
* `Contact`      (cycle parameterizing contact curves)
* `R1`           (first derived functor of direct image of ``ev^*O(-k)``)
"""
function AtiyahBottFormula(n::Int64, deg::Int64, n_marks::Int64, P, do_check::Bool = true, show_bar::Bool = true, down_col::Bool = true)::Vector{Rational{BigInt}}
    
    if n < 1
        printstyled("ERROR: ", bold=true, color=:red)
        println("n must be positive, correct ", n)
        return [Rational{BigInt}(0)]
    end
    if deg > 13 || deg < 1
        printstyled("ERROR: ", bold=true, color=:red)
        println("d must be between 1 and 13, correct ", deg)
        return [Rational{BigInt}(0)]
    end
    if n_marks < 0
        printstyled("ERROR: ", bold=true, color=:red)
        println("m must be non negative, correct ", n_marks)
        return [Rational{BigInt}(0)]
    end
    if !isa(P, Array)  #we want that P is an array, possibly with only one element
        P = [P]
    end
    if do_check
        if !is_zero_cycle(n, deg, n_marks, P)
            return [Rational{BigInt}(0)]
        end
    end

    if down_col && !fill_Data(n, deg)
        printstyled("ERROR: ", bold=true, color=:red)
        println("I was unable to download the colorations from the repository. Check your internet connexion or execute:")
        printstyled("julia> ", bold=true, color=:light_green)
        println("AtiyahBottFormula($n,$deg,$n_marks,P,$do_check,$show_bar,false);")
        return [Rational{BigInt}(0)]
    end

    
    local n_results::Int64 = length(P)   #this store the number of final results of our computation
    local result::Vector{Rational{BigInt}} = [Rational{BigInt}(0) for _ in 1:n_results] #the array of final results
    local max_col::Int64 = n+1   #the colors are number from 1 to n+1
     
    list_g::IOStream = open(current_dir*"/list_trees.txt") 
    #open the file containing the list of Prufer sequences of graphs

    s::Vector{Rational{BigInt}} = convert(Vector{Rational{BigInt}},rand(-1000*max_col*deg:1000*max_col*deg,max_col)) 

    #set up progress data                                    
    threshold::Int64 = sum(v -> number_trees[v-1]*max_col*(n^(v-1))*(v^n_marks), 2:deg+1)


    progress_data = ProgressData( 
        Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green, enabled=show_bar), #progress_bar
        0, #current_graph
        0, #top_aut
        threshold) #threshold

    for v in 2:(deg+1) #run the computation among all graphs with fixed number of vertices
                
        n_trees_nv = number_trees[v - 1]  #we known how many graphs there are with fixed number of vertices
        
        for n_g in 1:n_trees_nv  #run the computation for a fixed graph

            str = readline(list_g) #read a new line, we expect a Prufer seq plus the number of automorphisms
            local (g, aut) = get_graph(str)  #g is the graph, aut is the number of automorphisms
            progress_data.top_aut = aut
            #local top_aut::Int64 = copy(aut) #the number of automorphisms of g without any decoration

            pruf_str = string(split(str, ',')[1])
            result += AtiyahBottFormulaForGraph(g, pruf_str, aut, n, deg, 
                    n_marks, P, s, 
                    progress_data)

        end
    end
    
    close(list_g)     #close the file with Prufer sequences
    if n_results == 1
        println("Result: ", result[1])
    else 
        for res in 1:n_results
            println("Result number ", res, ": ", result[res])
        end
    end
    return result
end
