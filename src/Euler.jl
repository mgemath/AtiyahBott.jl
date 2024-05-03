# """
#     Euler_inv(g, c, w, s, m)

# The inverse of the (equivariant) Euler class of the normal bundle. This function is invoked automatically.
# # Arguments
# - `g::SimpleGraph{Int64}`: the graph.
# - `c::Vector{UInt8}`: the coloration.
# - `w::Vector{Int64}`: the weights.
# - `s::Tuple{Vararg{fmpq}}`: the scalars.
# - `m::Vector{Int64}`: the marks.

# """
function Euler_inv(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, mark::Marks_type, otd::Dict{Int64, fmpq})::fmpq
   
    local V::fmpq = fmpq(1)
    #local E::fmpq = fmpq(1)
    local temp1::fmpq = fmpq(1)
    # local temp2::fmpq = fmpq(1)
    local q1::fmpq
    # local p1::fmpq    
    local s1::fmpq  
    # col = Dict(vertices(g).=> col) #assing colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    #=
    omega_inv = Dict(edges(g).=> [d[e]//(scalars[col[src(e)]]-scalars[col[dst(e)]]) for e in edges(g)]) 
    merge!(omega_inv,Dict(reverse.(edges(g)).=> [d[e]//(scalars[col[dst(e)]]-scalars[col[src(e)]]) for e in edges(g)]))
    =#
    max_col = length(scalars)
    
    # for e in edges(g)
    #     q1 = fmpq(1)
    #     for j in 1:max_col
    #         if j != col[src(e)] && j != col[dst(e)]
    #             for alph in 0:d[e]
    #                 eq!(temp1, scalars[col[src(e)]])
    #                 eq!(temp2, scalars[col[dst(e)]])
    #                 mul_eq!(temp1, alph)
    #                 mul_eq!(temp2, d[e]-alph)
    #                 add_eq!(temp1, temp2)
    #                 div_eq!(temp1, d[e])
    #                 sub!(temp1, scalars[j])
    #                 mul_eq!(q1, temp1)
    #                 #q1 *= ((alph*scalars[col[src(e)]]+(d[e]-alph)*scalars[col[dst(e)]])//d[e]-scalars[j])
    #             end
    #         end
    #     end

    #     #eq!(temp1, omega_inv[e])
    #     eq!(temp1, scalars[col[dst(e)]])
    #     neg!(temp1)
    #     add_eq!(temp1, scalars[col[src(e)]])
    #     inv!(temp1)
    #     mul_eq!(temp1, d[e])

    #     pow_eq!(temp1, 2*d[e])
    #     if isodd(d[e])
    #         neg!(temp1)
    #     end
    #     q1 *= Nemo.factorial(fmpz(d[e]))^2
    #     # mul_eq!(q1, Nemo.factorial(fmpz(d[e]))^2)
    #     #div_eq!(temp1, factorial(d[e])^2)
    #     div_eq!(temp1, q1)
    #     mul_eq!(V, temp1)
    #     #E *= ((omega_inv[e])^(2*d[e]))*((-1)^d[e])//(factorial(d[e])^2)//q1
    # end
    
    for v in vertices(g)
        nghbrs = all_neighbors(g, v)
        # p1 = fmpq(1)
        # for j in 1:max_col
        #     if j != col[v]
        #         eq!(temp1, scalars[col[v]])
        #         sub!(temp1, scalars[j])
        #         mul_eq!(p1, temp1)
        #         #p1 *= scalars[col[v]]-scalars[j]
        #     end
        # end
        # pow_eq!(p1[col[v]], length(nghbrs)-1)
        # mul_eq!(V, p1)
        #p1 ^= length(nghbrs)-1
        mul_eq!(V, otd[col[v]]^(length(nghbrs)-1))
        
        s1 = fmpq(0)
        #local r1 = fmpq(1)
        
        for w in nghbrs
            eq!(temp1, scalars[col[w]])
            neg!(temp1)
            add_eq!(temp1, scalars[col[v]])
            inv!(temp1)
            # e = SimpleEdge(v,w)
            # d_e = haskey(d,e) ? d[e] : d[reverse(e)]
            # mul_eq!(temp1, d_e)
            mul_eq!(temp1, d[SimpleEdge(min(v,w),max(v,w))])
            # e = SimpleEdge(v,w)
            add_eq!(s1, temp1)
            mul_eq!(V, temp1)
            # s1 += omega_inv[e]
            # r1 *= omega_inv[e]
        end
        
        pow_eq!(s1, length(nghbrs) + num_marks(mark,v) - 3)
        #s1 ^= length(nghbrs) + num_marks(mark,v) - 3
        #mul_eq!(V, p1)
        mul_eq!(V, s1)
        #V *= p1*s1*r1

    end
    
    #mul_eq!(V, E)
    return V
end

function Lambda_Gamma_e(scalars::Tuple{Vararg{fmpq}}, d_e::Int64, col_1::Int64, col_2::Int64)::fmpq

    local V::fmpq = fmpq(1)
    #local E::fmpq = fmpq(1)
    local temp1::fmpq = fmpq(1)
    local temp2::fmpq = fmpq(1)
    max_col = length(scalars)

    q1 = fmpq(1)
    for j in 1:max_col
        if j != col_1 && j != col_2
            for alph in 0:d_e
                eq!(temp1, scalars[col_1])
                eq!(temp2, scalars[col_2])
                mul_eq!(temp1, alph)
                mul_eq!(temp2, d_e-alph)
                add_eq!(temp1, temp2)
                div_eq!(temp1, d_e)
                sub!(temp1, scalars[j])
                mul_eq!(q1, temp1)
                #q1 *= ((alph*scalars[col_1]+(d_e-alph)*scalars[col_2])//d_e-scalars[j])
            end
        end
    end

    #eq!(temp1, omega_inv[e])
    eq!(temp1, scalars[col_2])
    neg!(temp1)
    add_eq!(temp1, scalars[col_1])
    inv!(temp1)
    mul_eq!(temp1, d_e)

    pow_eq!(temp1, 2*d_e)
    if isodd(d_e)
        neg!(temp1)
    end
    q1 *= Nemo.factorial(fmpz(d_e))^2
    # mul_eq!(q1, Nemo.factorial(fmpz(d_e))^2)
    #div_eq!(temp1, factorial(d_e)^2)
    div_eq!(temp1, q1)
    mul_eq!(V, temp1)
    #E *= ((omega_inv[e])^(2*d_e))*((-1)^d_e)//(factorial(d_e)^2)//q1
    return V
end