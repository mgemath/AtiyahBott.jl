export
    AtiyahBottFormula


"""
    AtiyahBottFormula(n, d, m, P; do_check, show_bar)

Apply the Atiyah-Bott residue formula to the class `P`, in the moduli space of rational marked stable maps to the projective space of dimension `n` of degree `d` with `m` marks.
# Arguments
- `n::Int64`: the dimension of the projective space.
- `d::Int64`: the degree of the stable maps, it must be positive.
- `m::Int64`: the number of marks.
- `P`: the equivariant class.
- `do_check::Bool`: if `true`, checks if `P` is a well defined zero cycle, and stops the computation if this is not true. If `false`, the computation may have an unexpected behaviour. By default is `true`.
- `show_bar::Bool`: hide the progress bar if and only if this condition is `false`. By default is `true`.

The general construction of `P` is the following:
```julia-repl
julia> P = 
```
After `=`, one has to write an expression in the equivariant classes. The expression is a combination of the equivariant classes. We compute the degree of `P` by

```julia-repl
julia> AtiyahBottFormula(n,d,m,P);
```

# Example
```julia-repl
julia> P = Hypersurface(5);
julia> AtiyahBottFormula(3,1,0,P);
Warning: the class is not a 0-cycle.
julia> AtiyahBottFormula(4,1,0,P);
Result: 2875
julia> AtiyahBottFormula(4,1,0,P,do_check = false); #skip the preliminary check on `P`
julia> AtiyahBottFormula(4,1,0,P,show_bar = false); #it does not show the progress bar
```

The function returns an array of the same dimension of `P` (non-vectorized classes are assumed as 1-dimensional arrays). The Julia notation for accessing to array is `name_of_array[i]` where `i` is an index starting from 1.

# Example
```julia-repl
julia> P = Incidency(2)*Hypersurface(3);
julia> x = AtiyahBottFormula(3,2,0,P)[1];
Result: 81
julia> x
81
```

The class `P` supports parameters.
```julia-repl
julia> P = Hypersurface(3)*(Incidency(2)//3)^(d-1);
julia> d = 2;
julia> AtiyahBottFormula(3,d,0,P);
Result: 27
julia> d = 3;
julia> AtiyahBottFormula(3,d,0,P);
Result: 84
```

More examples are available in the support of the equivariant classes. It is enough to type `?` and then the name of the class. Currently, the supported classes are:

* `O1_i(j)`         (Euler class of ``\\mathrm{ev}_j^*\\mathcal{O}_{\\mathbb{P}^n}(1)``)
* `O1()`            (product of all `O1_i`)
* `Psi(a)`          (cycle of ``\\psi``-classes)
* `Jet(p,q)`        (Euler class of the jet bundle ``J^p`` of ``\\mathrm{ev}^*\\mathcal{O}_{\\mathbb{P}^n}(q)``)
* `Hypersurface(b)` (Euler class of the direct image of ``\\mathrm{ev}^*\\mathcal{O}_{\\mathbb{P}^n}(b)``)
* `Incidency(r)`    (cycle parameterizing curves meeting a linear subspace of codimension ``r``)
* `Contact()`       (cycle parameterizing contact curves)
* `R1(k)`           (first derived functor of direct image of ``\\mathrm{ev}^*\\mathcal{O}_{\\mathbb{P}^n}(-k)``)

To add more classes, please contact the authors.
"""
function AtiyahBottFormula(n::Int64, deg::Int64, n_marks::Int64, P_input; do_check::Bool = true, show_bar::Bool = true)::Vector{fmpq}
    
    if n < 1
        printstyled("ERROR: ", bold=true, color=:red)
        println("n must be positive, correct ", n)
        return [fmpq(0)]
    end
    if deg < 1 # deg > 13 || deg < 1
        printstyled("ERROR: ", bold=true, color=:red)
        println("d must be positive, correct ", deg)
        return [fmpq(0)]
    end
    if n_marks < 0
        printstyled("ERROR: ", bold=true, color=:red)
        println("m must be non negative, correct ", n_marks)
        return [fmpq(0)]
    end
    
    local n_results::Int64 = 1

    if isa(P_input, Array)
        n_results = length(P_input)
    end

    local P::Vector{Function} = Vector(undef, n_results)
    
    if isa(P_input, Array)
        for i in eachindex(P)
            P[i] = P_input[i].func
        end
    else
        P[1] = P_input.func
    end
    
    if do_check && !is_zero_cycle(n, deg, n_marks, P)
        return [fmpq(0)]
    end
    
    local result::Vector{Vector{fmpq}} = [[fmpq() for _ in 1:n_results] for _ in 1:Threads.nthreads()]
    local s::NTuple{n+1, fmpq} = (fmpq.(rand(Int16, n+1))...,)
    

    if show_bar #set up progress data
        number_trees = A000055(deg+1)
        threshold = sum(v -> number_trees[v]*(n+1)*(n^(v-1))*(v^n_marks), 2:deg+1)
        progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
        current_graph::Threads.Atomic{Int64} = Threads.Atomic{Int}(0)
    end
    
    last_ne::Int64 = 1
    all_weights::Vector{Vector{Int64}} = [[deg]]
    
    for ls in Iterators.flatten([TreeIt(v) for v in 2:(deg+1)])
        
        g::SimpleGraph{Int64} = LStoGraph(ls)

        if ne(g) > last_ne
            all_weights = get_weights(ne(g), deg)
            last_ne = ne(g)
        end


        top_aut::Int64 = count_iso(ls)

        l = Threads.SpinLock()
    
        Threads.@threads for c in collect(ColorsIt(ls, n+1)) #cols   #we run among all colorations of g

            local aut::Int64 = count_iso(ls, c)

            for m in Base.Iterators.product(repeat([1:nv(g)], n_marks)...)    #we run among all marks of g, if n_marks==0 we have only the empty mark                           
                
                for w in all_weights #we run among all weights of g
                    
                    try
                        local Euler::fmpq = fmpq(0)
                        
                        eq!(Euler, Euler_inv(g,c,w,s,m))
                        div_eq!(Euler, aut*prod(w))
                                                    
                        for res in 1:n_results      #compute each term of the array P
                            local temp::fmpq = fmpq(0)
                            eq!(temp, Base.invokelatest(P[res], g,c,w,s,m))
                            # eq!(temp, P[res](g,c,w,s,m))
                            mul_eq!(temp, Euler)
                            add_eq!(result[Threads.threadid()][res], temp)
                            # result[res] += P[res](g,c,w,s,m)*Euler    #apply Atiyah-Bott
                        end
                        
                    catch err 
                        if isa(err, DivideError) 
                            error("Some division by zero occurred. Try again")
                        end
                        println(err)
                        error("Some error occurred")
                        return zeros(fmpq, n_results)
                    end
                end
                
                if show_bar
                    Threads.atomic_add!(current_graph, top_aut÷aut)
                    #progress_data.current_graph += progress_data.top_aut÷aut   
                    #update the progress bar
                    Threads.lock(l)
                    update!(progress_bar, current_graph[],
                            showvalues = [(:"Total number of graphs",threshold),(:"Current graph",current_graph[])])
                    Threads.unlock(l)
                end
            end
        end
    end
    
    for nt in 2:Threads.nthreads()
        add_eq!(result[1], result[nt])
    end
    
    if n_results == 1
        println("Result: ", result[1][1])
    else 
        for res in 1:n_results
            println("Result number ", res, ": ", result[1][res])
        end
    end
    
    return result[1]
end
