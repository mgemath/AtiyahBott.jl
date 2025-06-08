function Euler_inv(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, mark::Marks_type, otd::Dict{Int64, fmpq})::fmpq
   
    local V::fmpq = one(scalars[1])
    d = Dict(Graphs.edges(g).=> weights) #assign weights to edges
    max_col = length(scalars)
    
    for v in Graphs.vertices(g)
        nghbrs = Graphs.all_neighbors(g, v)
        V *= otd[col[v]]^(length(nghbrs)-1)
        
        s1 = zero(scalars[1])
        
        for w in nghbrs
            temp1 = d[SimpleEdge(min(v,w),max(v,w))]//(scalars[col[v]] - scalars[col[w]])
            s1 += temp1
            V *= temp1
        end
        s1 ^= length(nghbrs) + num_marks(mark,v) - 3
        V *= s1
    end
    return V
end

function Lambda_Gamma_e(scalars::Tuple{Vararg{fmpq}}, d_e::Int64, col_1::Int64, col_2::Int64)::fmpq

    max_col = length(scalars)

    q1 = one(scalars[1])
    for j in 1:max_col
        if j != col_1 && j != col_2
            for alph in 0:d_e
                q1 *= d_e//(alph*scalars[col_1]+(d_e-alph)*scalars[col_2]-d_e*scalars[j])
            end
        end
    end

    return q1* (((-1)^d_e) // ((Nemo.factorial(fmpz(d_e)))^2) )*( (d_e//(scalars[col_1] - scalars[col_2]))^(2*d_e) )
end