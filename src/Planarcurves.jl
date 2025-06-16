export hypersur, planarcurves, computeplanarcurves, test_P_n_dual

function planarcurves(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, s::Tuple{Vararg{fmpq}}, r::Tuple{Vararg{fmpq}}, fixed_point::Int64)::fmpq

    local p1::fmpq = one(s[1])
    d = Dict(Graphs.edges(g).=> weights) #assign weights to edges
    
    for e in Graphs.edges(g)
        for alph in 0:(d[e])
            p1 *= (alph*s[col[Graphs.src(e)]]+(d[e]-alph)*s[col[Graphs.dst(e)]] + r[fixed_point])//d[e]
        end
    end
    
    for v in Graphs.vertices(g)
        p1 *= (s[col[v]] + r[fixed_point])^(1-length(Graphs.all_neighbors(g, v)))   
    end

    return p1

end

function planarcurves( b , r, f)::EquivariantClass
    
    rule = :(planarcurves( g, c, w, s, $b, $r, $f ))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end

function computeplanarcurves(n::Int64, deg::Int64, n_marks::Int64, b::Int64; show_bar::Bool = true)::Vector{fmpq}

    R, _s, _r = polynomial_ring(QQ, :s => 1:(n+1), :r => 1:(n+1))
    S = fraction_field(R)
    s = ([S(_s[i]) for i in 1:(n+1)]...,)  # lambda
    r = ([S(_r[i]) for i in 1:(n+1)]...,)  # mu
    
    # if n < 1
    #     printstyled("ERROR: ", bold=true, color=:red)
    #     println("n must be positive, correct ", n)
    #     return [zero(s[1])]
    # end
    # if deg < 1 # deg > 13 || deg < 1
    #     printstyled("ERROR: ", bold=true, color=:red)
    #     println("d must be positive, correct ", deg)
    #     return [zero(s[1])]
    # end
    # if n_marks < 0
    #     printstyled("ERROR: ", bold=true, color=:red)
    #     println("m must be non negative, correct ", n_marks)
    #     return [zero(s[1])]
    # end
    
    local n_results::Int64 = 1

    # if isa(P_input, Array)
    #     n_results = length(P_input)
    # end

    # local P::Vector{Function} = Vector(undef, n_results)
    
    # if isa(P_input, Array)
    #     for i in eachindex(P)
    #         P[i] = P_input[i].func
    #     end
    # else
    #     P[1] = P_input.func
    # end
    
    local result::Vector{Vector{fmpq}} = [[zero(s[1]) for _ in 1:n_results] for _ in 1:Threads.nthreads()]

    nc = Dict{Int64,Vector{Int64}}([i for i in 1:(n+1)] .=> [[j + Int64(i<=j) for j in 1:n] for i in 1:(n+1)])
    Lambda_Gamma_e_dict::Dict{Tuple{Int64, Int64, Int64}, fmpq} = Dict{Tuple{Int64, Int64, Int64}, fmpq}()
    omega_t_dict::Dict{Int64, fmpq} = Dict{Int64, fmpq}()
    for c_1 in 1:(n+1)
        omega_t_dict[c_1] = one(s[1])
        for c_2 in 1:(n+1)
            if c_2 > c_1
                for deg_e in 1:deg
                    Lambda_Gamma_e_dict[deg_e, c_1, c_2] = Lambda_Gamma_e(s, deg_e, c_1, c_2)
                end
            end
            if c_2 != c_1
                omega_t_dict[c_1] *= s[c_1] - s[c_2]
            end
        end
    end
    

    if show_bar #set up progress data
        number_trees = A000055(deg+1)
        threshold = sum(v -> number_trees[v]*(n+1)*(n^(v-1))*binomial(v+n_marks-1,n_marks), 2:deg+1)
        progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
        current_graph::Threads.Atomic{Int64} = Threads.Atomic{Int}(0)
    end
    
    last_ne::Int64 = 1
    all_weights::Vector{Vector{Int64}} = [[deg]]
    
    for ls in Iterators.flatten([TreeIt(v) for v in 2:(deg+1)])
        
        g::SimpleGraph{Int64} = LStoGraph(ls)

        if Graphs.ne(g) > last_ne
            all_weights = get_weights(Graphs.ne(g), deg)
            last_ne = Graphs.ne(g)
        end


        tree_aut::Int64 = count_iso(ls)

        CI, parents, subgraph_ends = col_it_init(ls, nc)
        for col in collect(CI)

            local top_aut::Int64 = count_iso(ls, col)

            for m_inv in with_replacement_combinations(1:Graphs.nv(g), n_marks)
            # for m in Base.Iterators.product(repeat([1:nv(g)], n_marks)...)    #we run among all marks of g, if n_marks==0 we have only the empty mark                           
                aut = count_iso(ls, col, m_inv)
                for w in all_weights #we run among all weights of g
                    PRODW = prod(w)
                    d = Dict(Graphs.edges(g).=> w)
                    try
                        local Euler::fmpq = zero(s[1])
                        local temp = Vector{fmpq}(undef, n_results)
                        
                        for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, n_marks))

                            # for res in eachindex(temp)
                            #     # temp[res] = Base.invokelatest(P[res], g, col, w, s, m) * hypersur(g, col, w, s, b)
                            #     temp[res] = hypersur(g, col, w, s, b)*(Incidency(g, col, w, r, n)^2)
                            # end

                            # all(res -> temp[res] == zero(s[1]), eachindex(temp)) && continue # check if at least one partial result is not zero
                            
                            if Euler == zero(s[1])
                                Euler = Euler_inv(g, col, w, s, m, omega_t_dict)//(aut*PRODW)

                                for e in Graphs.edges(g)
                                    triple = (d[e], min(col[Graphs.src(e)], col[Graphs.dst(e)]), max(col[Graphs.src(e)], col[Graphs.dst(e)]))

                                    Euler *= Lambda_Gamma_e_dict[triple]
                                end
                            end

                            for fixed_p in 1:(n+1)
                                # temp[1] = (hypersur(g, col, w, s, b)  * Euler) # M_bar part
                                # temp[1] *= (r[fixed_p]^n)*(Euler_2(r, fixed_p)) # P^n dual part

                                temp[1] =  r[fixed_p]*planarcurves(g, col, w, s, r, fixed_p) * Euler_2(r, fixed_p)
                                temp[1] *= Incidency(g, col, w, s, n)^2 * Euler
                                # computeplanarcurves(3, 1, 0, 0);

                                # temp[1] =  planarcurves(g, col, w, s, r, fixed_p) * Euler_2(r, fixed_p)
                                # temp[1] *= Incidency(g, col, w, s, n)^3 *Incidency(g, col, w, s, n-1)^2 * Euler

                                result[1][1] += temp[1]
                            end

                            # temp[1] = (hypersur(g, col, w, s, b)  * Euler) # M_bar part
                            # result[1][1] += temp[1]

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
                    Threads.atomic_add!(current_graph, tree_aut÷top_aut)
                    #progress_data.current_graph += progress_data.tree_aut÷top_aut   
                    #update the progress bar
                    update!(progress_bar, current_graph[],
                            showvalues = [(:"Total number of graphs",threshold),(:"Current graph",current_graph[])])
                end
            end
        end
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

function Euler_2(r, fixed_p::Int64)

    p1 = one(r[1])

    for j in 1:length(r)
        j == fixed_p && continue
        p1 *= 1//(r[fixed_p] - r[j])
    end

    return p1    
end

function test_P_n_dual(n::Int64)

    R, _s = polynomial_ring(QQ, :x => 1:((n+1)))
    S = fraction_field(R)
    r = ([S(_s[i]) for i in 1:(n+1)]...,)

    p1 = zero(r[1])

    for fixed_p in 1:length(r)
        p1 += (r[fixed_p]^n)*(Euler_2(r, fixed_p))
    end

    return p1
    
end

function hypersur(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, s::Tuple{Vararg{fmpq}}, b::Int64)::fmpq

    local p1::fmpq = one(s[1])
    d = Dict(Graphs.edges(g).=> weights) #assign weights to edges
    
    for e in Graphs.edges(g)
        for alph in 0:(b*d[e])
            p1 *= (alph*s[col[Graphs.src(e)]]+(b*d[e]-alph)*s[col[Graphs.dst(e)]])//d[e]
        end
    end
    
    for v in Graphs.vertices(g)
        p1 *= (b*s[col[v]])^(1-length(Graphs.all_neighbors(g, v)))   
    end

    return p1

end

function hypersur( b )::EquivariantClass
    
    rule = :(hypersur( g, c, w, s, $b ))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end