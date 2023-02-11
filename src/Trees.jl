# This is an implementation of the algorithm to generate all trees on n vertices up to 
# isomorphism. The algorithm is described in the paper 
# CONSTANT TIME GENERATION OF FREE TREES*
# ROBERT ALAN WRIGHTS’, BRUCE RICHMONDT, ANDREW ODLYZKO AND BRENDAN D. MCKAY
#
# written by Csaba Schneider, contribution from Giosuè Muratore
struct TreeIt
    n_vert::Int64
end

function Base.eltype(::Type{TreeIt})
    Vector{Int64}
end

# in order to have an iterator that produces SimpleGraph, use: Base.Generator(LStoGraph, TreeIt(2))

Base.@propagate_inbounds function Base.length(TI::TreeIt)::Int64

    TI.n_vert < 4 && return 1
    a =  A000081(TI.n_vert) 
    
    len = a[TI.n_vert] - sum(i -> a[i]*a[TI.n_vert-i], 1:(TI.n_vert >>> 1))

    if iseven(TI.n_vert)
        a_n = a[TI.n_vert >>> 1]
        len += ((a_n + 1)*a_n) >>> 1
    end

    return len
end

function Base.iterate( TI::TreeIt, state::Tuple = () )::Union{Nothing, Tuple{Vector{Int64}, Tuple}}

    if isempty(state) #first iteration
        
        if TI.n_vert < 4
            return map(i -> min(2,i), 1:TI.n_vert), (0,)
            # return StarGraph(TI.n_vert), Dict('q' => 0)
        end

        n = TI.n_vert
        r = n÷2 + 1
        l = [ 1:r; 2:n-r+1 ];
        w = [ 0:r-1; 1:n-1]
        p, q, h1, h2 =  n, n-1, r, n
        
        if isodd( n ) 
            c = Inf
        else 
            c = Float64(n+1)
        end 

        return l, (q, l, w, n, p, h1, h2, c, r, false, false, false, false)
    end

    if state[1] == 0 #last iteration
        return nothing
    end

    # next
    (q, l, w, n, p, h1, h2, c, r, fixit, needr, needc, needh2) = state


    # fixit = false

    if c == n+1 || p == h2 && (l[h1] == l[h2]+1 && n-h2>r-h1 ||
            l[h1] == l[h2] && n-h2+1 < r - h1)  
        if l[r] > 3 
            p = r;  q = w[r]
            if h1 == r 
                h1 = h1 - 1
            end 
            fixit = true
        else 
            p = r; r = r-1; q = 2
        end
    end
    
    # needr = false; needc = false; needh2 = false
    if p <= h1 h1 = p - 1; end
    if p <= r 
        needr = true
    elseif p <= h2 
        needh2 = true 
    elseif l[h2] == l[h1] - 1 && n - h2 == r - h1 
        if p <= c needc = true end
    else 
        c = Inf
    end

    oldp = p; δ = q - p; oldlq = l[q]; oldwq = w[q]; p = Inf
    
    for i in oldp:n
        l[i] = l[i+δ]
        if l[i] == 2 
            w[i] = 1
        else
            p = i
            if l[i] == oldlq 
                q = oldwq
            else 
                q = w[i+δ] - δ
            end
            w[i] = q
        end
    
        if needr && l[i] == 2 
            needr = false; needh2 = true; r = i - 1
        end
        
        if needh2 && l[i] <= l[i-1] && i > r + 1 
            needh2 = false; h2 = i - 1
            if l[h2] == l[h1] - 1 && n - h2 == r - h1 
                needc = true
            else 
                c = Inf
            end
        end
        
        if needc 
            if l[i] != l[h1-h2+i] - 1
                needc = false; c = Float64(i)
            else 
                c = Float64(i+1)
            end
        end
    end

    if fixit 
        r = n-h1+1
        for i in (r+1):n 
            l[i] = i-r+1; w[i] = i-1
        end
        w[r+1] = 1; h2 = n; p = n; q = p-1; c = Inf
    else
        if p == Inf 
            if l[oldp-1] != 2 
                p = oldp - 1
            else 
                p = oldp - 2
            end
            q = w[p]
        end

        if needh2 
            h2 = n
            if l[h2] == l[h1] - 1 && h1 == r 
                c = Float64(n + 1)
            else 
                c = Inf
            end
        end
    end
    return l, (q, l, w, n, p, h1, h2, c, r, false, false, false, false)
end

function LStoGraph(ls::Vector{Int64})::SimpleGraph{Int64}

    n::Int64 = length( ls )
    ans::SimpleGraph{Int64} = SimpleGraph(n)

    for v in 2:n
        p = findlast(i -> i < v && ls[i] == ls[v] - 1, eachindex(ls))
        add_edge!(ans, p, v)
    end
    
    return ans
end

function A000081(n::Int64)::Vector{Int64} #OEIS A000081
    # a(n) = (1/(n-1)) * Sum_{k=1..n-1} ( Sum_{d|k} d*a(d) ) * a(n-k)
    n < 3 && return ones(Int64, n)

    a = Vector{Int64}(undef, n) #OEIS A000081
    a[1] = 1; a[2] = 1; a[3] = 2; a[4] = 4;

    local first::Int64
    local second::Int64

    for n in 5:n
        first = a[n-1]
        for k in 2:(n-1)
            second = k*a[k] + 1
            for d in 2:(k-1)
                k % d > 0 && continue
                second += d*a[d]
            end
            first += second*a[n-k]
        end
        a[n] = first÷(n-1)
    end

    return a
end

function A000055(n_vert::Int64)::Vector{Int64} # OEIS A000055

    n_vert < 4 && return ones(Int64, n_vert)
    
    a =  A000081(n_vert) #OEIS A000081
    len = Vector{Int64}(undef, n_vert) #OEIS A000055
    
    len[1] = 1
    foreach(j -> len[j]=a[j] - sum(i->a[i]*a[j-i], 1:(j >>> 1)), 2:n_vert)

    for j in 2:2:n_vert
        a_n = a[j >>> 1]
        len[j] += ((a_n + 1)*a_n) >>> 1
    end

    return len
end

# println(raw"This is an implementation of the algorithm to generate
# all trees on n vertices up to isomorphism.
# The algorithm is described in the paper:
#     CONSTANT TIME GENERATION OF FREE TREES
# by ROBERT ALAN WRIGHTS’, BRUCE RICHMONDT, ANDREW ODLYZKO AND BRENDAN D. MCKAY 
# written by Csaba Schneider, contribution from Giosuè Muratore.");
