export 
    Incidency,
    Hypersurface,
    Contact,
    O1_i,
    O1,
    R1,
    Psi,
    Jet

using LightGraphs: SimpleEdge

"""
Equivariant class of the cycle parameterizing curves meeting a linear subspace of codimension `r`.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `r::Int64`: the codimension of the subvariety. Alternatively, it can be an array of integers, meaning the multiplication of the equivariant class defined by each element of the array.

# Example
```julia-repl
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,3)^2;
julia> AtiyahBottFormula(3,1,0,P);
Result: 1//1
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,[2,2,3]);
julia> AtiyahBottFormula(3,1,0,P);
Result: 1//1
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,[2,2])*Hypersurface(g,c,w,s,3);
julia> AtiyahBottFormula(3,3,0,P);
Result: 756//1
```
!!! warning "Attention!"

    The program will stop if `r` is not positive.

"""
function Incidency(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, r::Int64)::Rational{BigInt}
    
    local p1 = Rational{BigInt}(0); #the final result
    r -= 1                          
    
    # col = Dict(vertices(g).=> coloration) #assign colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges

    for e in edges(g)
        for t in (0:r)
            p1 += d[e]*(scalars[col[src(e)]]^(t))*(scalars[col[dst(e)]]^(r-t))
        end
    end

    return p1

end

function Incidency(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, r::Vector{Int64})::Rational{BigInt}
    
    local p1 = Rational{BigInt}(1);
    
    for j in unique(r)
        p1 *= Incidency(g, coloration, weights, scalars, j)^count(x->x==j, r)
    end
    return p1

end


"""
Equivariant class of the Euler class of the bundle equal to the direct image under the forgetful map of ``ev^*O(b)``. It parameterizes curves contained in a hypersurface of degree `b`.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `b::Int64`: the degrees of the hypersurface. Alternatively, it can be an array of integers, meaning the multiplication of the equivariant class defined by each element of the array.

# Example
```julia-repl
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);
julia> AtiyahBottFormula(4,1,0,P);
Result: 2875//1
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,[3,3]);
julia> AtiyahBottFormula(5,2,0,P);
Result: 423549//8
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,4)*Hypersurface(g,c,w,s,2);
julia> AtiyahBottFormula(5,3,0,P);
Result: 422690816//27
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,2)^4;
julia> AtiyahBottFormula(7,4,0,P);
Result: 25705160//1
```
!!! warning "Attention!"

    The program will stop if `b` is not positive.

"""
function Hypersurface(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, b::Int64)::Rational{BigInt}

    local p1 = Rational{BigInt}(1)
    local q1 = Rational{BigInt}(1)
    
    # col = Dict(vertices(g).=> coloration) #assign colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    for e in edges(g)
        for alph in 0:(b*d[e])
            p1 *= (alph*scalars[col[src(e)]]+(b*d[e]-alph)*scalars[col[dst(e)]])//d[e]
        end
    end
    
    for v in vertices(g)
        q1 *= (b*scalars[col[v]])^(1-length(all_neighbors(g, v)))   
    end

    return p1*q1

end

function Hypersurface(g::SimpleGraph, coloration::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, b::Vector{Int64})::Rational{BigInt}

    local p1 = Rational{BigInt}(1)
    
    for j in unique(b)
        p1 *= Hypersurface(g, coloration, weights, scalars, j)^count(x->x==j, b)
    end
    return p1
end


"""
Equivariant class of the Euler class of the bundle equal to the direct image under the forgetful map of: ``ev^*O(2)`` tensor the dualizing sheaf of the forgetful map. It parameterizes contact curves in an odd dimensional projective space.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.

# Example
```julia-repl
julia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)^2*O1_i(g,c,w,s,m,2)^3*Contact(g,c,w,s);
julia> AtiyahBottFormula(3,1,2,P);
Result: 1//1
```
"""
function Contact(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}})::Rational{BigInt}

    local p1 = Rational{BigInt}(1)
    local q1 = Rational{BigInt}(1)
    
    # col = Dict(vertices(g).=> coloration) #assign colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    for e in edges(g)
        for alph in (1:2*d[e]-1)
            p1 *= (alph*scalars[col[src(e)]]+(2*d[e]-alph)*scalars[col[dst(e)]])//d[e]
        end
    end
    
    for v in vertices(g)
        q1 *= (2*scalars[col[v]])^(length(all_neighbors(g, v))-1)
    end

    return p1*q1

end


"""
Equivariant class of the pull-back of ``O(1)`` with respect to the i-th evaluation map.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `m::marks`: the marks.
- `i::Int64`: the evaluation map.

# Example
The following Gromov-Witten invariants
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^{2},1)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)\\cdot\\mathrm{ev}_{2}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1) &= 1 \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{3},1)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2}\\cdot\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(2))) &= 4
\\end{aligned}
```
can be computed as
```julia-repl
julia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)*O1_i(g,c,w,s,m,2);
julia> AtiyahBottFormula(2,1,2,P);
Result: 1//1
julia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)^2*Hypersurface(g,c,w,s,2);
julia> AtiyahBottFormula(3,1,1,P);
Result: 4//1
```
!!! warning "Attention!"

    The program will stop if `i` is not between 1 and the number of marks.

"""
function O1_i(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, mark::marks, i::Int64)::Rational{BigInt}
    
    # col = Dict(vertices(g).=> coloration) #assing colors to vertices
    
    return scalars[col[mark.get_vertex[i]]]
end


"""
Equivariant class of the pull-back of ``O(1)`` with respect to the product of all evaluation maps.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `m::marks`: the marks.

This function is equivalent to the product of the function `O1_i(g,c,w,s,m,i)` where `i` runs from 1 to the number of marks.

# Example
```julia-repl
julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m);
julia> AtiyahBottFormula(2,3,8,P);
Result: 12//1
julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Hypersurface(g,c,w,s,3);
julia> AtiyahBottFormula(3,2,1,P);
Result: 81//1
```

In order to remove `O1_i(g,c,w,s,m,i)` for some `i`, it is enough to divide by that function.

# Example
```julia-repl
julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)//O1_i(g,c,w,s,m,1);
```
Here `P` is the product of all `O1_i(g,c,w,s,m,i)` where `i` runs from 2 to `m`.
"""
function O1(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, mark::marks)::Rational{BigInt}
    
    local p1 = Rational{BigInt}(1)
    # col = Dict(vertices(g).=> coloration)
    
    for t in 1:mark.m
        p1 *= scalars[col[mark.get_vertex[t]]]
    end
    
    return p1
end

"""
The equivariant class of the first derived functor of the pull-back of ``O(-k)``.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `k::Int64`: a positive integer.

# Example
```julia-repl
julia> P = (g,c,w,s,m) -> R1(g,c,w,s,1)^2;
julia> AtiyahBottFormula(1,3,0,P);
Result: 1//27
```
!!! warning "Attention!"

    The program will stop if `k` is not positive.

"""
function R1(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, k::Int64)::Rational{BigInt}
    
    local p1 = Rational{BigInt}(1)
    local q1 = Rational{BigInt}(1)
    
    # col = Dict(vertices(g).=> coloration) #assign colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    for e in edges(g)
        for alph in 1:(k*d[e]-1)
            p1 *= (-1)*(alph*scalars[col[src(e)]]+(k*d[e]-alph)*scalars[col[dst(e)]])//d[e]
        end
    end
    
    for v in vertices(g)
        q1 *= (-k*scalars[col[v]])^(length(all_neighbors(g, v))-1)   
    end

    return p1*q1

end

"""
Equivariant class of the cycle of ``psi``-classes.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `m::marks`: the marks.
- `a::Vector{Int64}`: the vector of the exponents of the ``psi`` classes. It is ordered, meaning that the first element is the exponent of ``psi_1``, the second is the exponent of ``psi_2``, and so on.

!!! note

    The size of `a` must be at most `m`. If it is smaller, missing exponents will be considered as zeros.
    If `a` is a number, it will be considered as the exponent of ``psi_1``.

!!! warning "Attention!"

    The program will stop if we have one of the following conditions:

    * the size of `a` is bigger than `m`,
    * `a` contains a negative number.

# Example
```julia-repl
julia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)^5*O1_i(g,c,w,s,m,2)^2*Hypersurface(g,c,w,s,5)*Psi(g,c,w,s,m,[1,0]);
julia> AtiyahBottFormula(6,2,2,P);
Result: 495000//1
julia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)^8*O1_i(g,c,w,s,m,2)^6*Hypersurface(g,c,w,s,7)*Psi(g,c,w,s,m,2);
julia> AtiyahBottFormula(10,2,2,P);
Result: 71804533752//1
julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Psi(g,c,w,s,m,4);
julia> AtiyahBottFormula(2,2,1,P);
Result: 1//8
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^4*O1_i(g,c,w,s,m,1)*(O1_i(g,c,w,s,m,1) + Psi(g,c,w,s,m,1))
julia> AtiyahBottFormula(2,2,1,P); #number of plane conics through four points and tangent to a line
Result: 2
```
!!! warning "Psi is singleton!"

    `Psi` cannot be multiplied by itself.
    ```julia-repl
    julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Psi(g,c,w,s,m,1)^4;                  #this is **wrong**
    julia> AtiyahBottFormula(2,2,1,P);
    Warning: more instances of Psi has been found. Type:
    julia> ?Psi
    for support.
    julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Psi(g,c,w,s,m,3)*Psi(g,c,w,s,m,1);   #this is **wrong**
    julia> AtiyahBottFormula(2,2,1,P);
    Warning: more instances of Psi has been found. Type:
    julia> ?Psi
    for support.
    julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Psi(g,c,w,s,m,4);
    julia> AtiyahBottFormula(2,2,1,P);
    Result: 1//8
    julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*(Psi(g,c,w,s,m,3)*O1(g,c,w,s,m)+Psi(g,c,w,s,m,2)*O1(g,c,w,s,m)^2);
    julia> AtiyahBottFormula(2,2,1,P);
    Result: 1//8
    julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)*(Psi(g,c,w,s,m,7)*O1(g,c,w,s,m)+Psi(g,c,w,s,m,6)*O1(g,c,w,s,m)^2);
    julia> AtiyahBottFormula(3,2,1,P);
    Result: -5//16
    ```
"""
function Psi(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, mark::marks, a::Int64)::Rational{BigInt}
    
    return Psi(g, col, weights, scalars, mark, [a])
end
function Psi(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, mark::marks, a::Vector{Int64})::Rational{BigInt}
    
    if findfirst(x -> x>0, a) === nothing #if all of them are zero
        return Rational{BigInt}(1)
    end
    
    local q1 = Rational{BigInt}(1)
    
    for v in vertices(g)
        S_v = num_marks(mark, v)
        S_v == 0 && continue
        
        nghbrs = all_neighbors(g, v)
        N = S_v + length(nghbrs)
        d = Dict(edges(g).=> weights) #assign weights to edges
        
        if N == 2
            i = invert_marks(mark)[v][1]
            (i > length(a) || a[i] == 0) && continue
            w = nghbrs[1]
            e = SimpleEdge(v,w)
            d_e = haskey(d,e) ? d[e] : d[reverse(e)]
            q1 *= ((scalars[col[w]]-scalars[col[v]])//d_e)^a[i]
            continue
        end
        #If no previous condition holds, then N>2
        if length(a) < mark.m #correct the dimension of a, if it is too low
            a = vcat(a,[0 for _ in 1:(mark.m-length(a))])
        end
        a_v = [a[i] for i in invert_marks(mark)[v]]
        Sum_ai = sum(a_v)
        if Sum_ai > N-3
            return Rational{BigInt}(0)
        end

        local s1 = Rational{BigInt}(0)

        omega_inv = Dict(edges(g).=> [d[e]//(scalars[col[src(e)]]-scalars[col[dst(e)]]) for e in edges(g)]) 
        merge!(omega_inv,Dict(reverse.(edges(g)).=> [d[e]//(scalars[col[dst(e)]]-scalars[col[src(e)]]) for e in edges(g)]))

        for w in nghbrs
            s1 += omega_inv[SimpleEdge(v,w)]
        end
        s1 ^= -Sum_ai
        q1 *= multinomial(N-3-Sum_ai,a_v...)*s1

        end
        
    return q1
end
"""
Equivariant class of the jet bundle ``J^p`` of the pull back of ``O(q)`` with respect to the first ``psi``-class.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `m::marks`: the marks.
- `p::Int64`: the exponent of the Jet bundle. In particular, it is a bundle of rank p+1.
- `q::Int64`: the degree of the line bundle that is pulled back.


!!! note

    In order to define this bundle, the number of marks must be at least 1.
    You cannot multiply this bundle by the class `Psi(g,c,w,s,m,a)`.


# Example
```julia-repl
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^4*Jet(g,c,w,s,m,1,1);
julia> AtiyahBottFormula(2,2,1,P);
Result: 2//1
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^4*(Jet(g,c,w,s,m,1,1)+O1(g,c,w,s,m)^2);
julia> AtiyahBottFormula(2,2,1,P);
Result: 3//1
julia> P = (g,c,w,s,m) -> (O1(g,c,w,s,m)^2)//k*Jet(g,c,w,s,m,4*d-2,k);
julia> d=1;k=1;AtiyahBottFormula(3,d,1,P);   #The value of this integral does not depend on k, only on d
```
"""
function Jet(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, mark::marks, p::Int64, q::Int64)::Rational{BigInt}
    
    local s1 = Rational{BigInt}(0)
    for h in 0:p
        s1 += stirlings1(p+1, p+1-h)*(q*O1_i(g, col, weights, scalars, mark,1))^(p+1-h)*Psi(g, col, weights, scalars, mark,[h])
    end

    return s1
end

"""
The inverse of the (equivariant) Euler class of the normal bundle.
# Arguments
- `g::SimpleGraph`: the graph.
- `c::Vector{UInt8}`: the coloration.
- `w::Vector{Int64}`: the weights.
- `s::Rational{BigInt}`: the scalars.
- `m::marks`: the marks.
"""
function Euler_inv(g::SimpleGraph, col::Vector{UInt8}, weights::Vector{Int64}, scalars::Vector{Rational{BigInt}}, mark::marks)::Rational{BigInt}
   
    local V = Rational{BigInt}(1)
    local E = Rational{BigInt}(1)
    
    # col = Dict(vertices(g).=> coloration) #assing colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    omega_inv = Dict(edges(g).=> [d[e]//(scalars[col[src(e)]]-scalars[col[dst(e)]]) for e in edges(g)]) 
    merge!(omega_inv,Dict(reverse.(edges(g)).=> [d[e]//(scalars[col[dst(e)]]-scalars[col[src(e)]]) for e in edges(g)]))
    
    max_col = length(scalars)
    
    for e in edges(g)
        local q1 = Rational{BigInt}(1)
        for j in 1:max_col
            if j != col[src(e)] && j != col[dst(e)]
                for alph in 0:d[e]
                    q1 *= ((alph*scalars[col[src(e)]]+(d[e]-alph)*scalars[col[dst(e)]])//d[e]-scalars[j])
                end
            end
        end
        E *= ((omega_inv[e])^(2*d[e]))*((-1)^d[e])//(factorial(d[e])^2)//q1
    end
    
    for v in vertices(g)
        nghbrs = all_neighbors(g, v)
        local p1 = Rational{BigInt}(1)
        for j in 1:max_col
            if j != col[v]
                p1 *= scalars[col[v]]-scalars[j]
            end
        end
        p1 = p1^(length(nghbrs)-1)
        
        local s1 = Rational{BigInt}(0)
        local r1 = Rational{BigInt}(1)
        
        for w in nghbrs
            e = SimpleEdge(v,w)
            s1 += omega_inv[e]
            r1 *= omega_inv[e]
        end
        s1 ^= length(nghbrs) + num_marks(mark,v) - 3
        V *= p1*s1*r1

    end

    return V*E
end

#############################
###List of the class useful for computing the rank
#############################
function Incidency(n::Int64, deg::Int64, c::Int64, s::Int64, r::Int64)::Cycle

    try
        r < 0 && error(string("Incidency requires a positive integer, correct ", r))
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end

    return Cycle(r-1, 0)
end

function Incidency(n::Int64, deg::Int64, c::Int64, s::Int64, r::Vector{Int64})::Cycle

    try
        if findfirst(x -> x<1, r) !== nothing #if some of them is non positive
            error(string("Incidency requires a positive integer, correct ", r))
        end
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end

    return Cycle(sum(r)-length(r), 0)
end


function Hypersurface(n::Int64, deg::Int64, c::Int64, s::Int64, b::Int64)::Cycle

    try
        b < 1 && error(string("Hypersurface requires a positive integer, correct ", b))
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    
    return Cycle(b*deg+1, 0)
end

function Hypersurface(n::Int64, deg::Int64, c::Int64, s::Int64, b::Vector{Int64})::Cycle

    try
        if findfirst(x -> x<1, b) !== nothing #if some of them is negative
            error(string("Hypersurface requires a positive integer, correct ", b))
        end
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end

    return Cycle(sum(b)*deg+length(b), 0)
end


function Contact(n::Int64, deg::Int64, c::Int64, s::Int64)::Cycle
    
    return Cycle(2*deg-1, 0)
end


function O1_i(n::Int64, deg::Int64, n_marks::Int64, c::Int64, s::Int64, i::Int64)::Cycle

    try
        if (i < 1 || i > n_marks)
            error(string("O1_i requires a positive integer between 1 and ", n_marks, ", correct ",i))
        end
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    
    return Cycle(1, 0)
end


function O1(n::Int64, deg::Int64, n_marks::Int64, c::Int64, s::Int64)::Cycle
    
    return Cycle(n_marks, 0)
end


function R1(n::Int64, deg::Int64, c::Int64, s::Int64, k::Int64)::Cycle

    try
        k < 1 && error(string("R1 requires a positive integer, correct ", k))
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end

    return Cycle(k*deg-1, 0)
end

function Psi(n::Int64, deg::Int64, n_marks::Int64, c::Int64, s::Int64, a::Int64)::Cycle
    
    try
        a < 0 && error(string("exponents of psi classes must be nonnegative, correct ", a))
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    psi_deg = a==0 ? 0 : 1

    return Cycle(a, psi_deg)
end

function Psi(n::Int64, deg::Int64, n_marks::Int64, c::Int64, s::Int64, a::Vector{Int64})::Cycle
    
    try
        length(a) > n_marks && error(string("size of ",a," is greater than ",n_marks))

        if findfirst(x -> x<0, a) !== nothing #if some of them is negative
            error(string("exponents of psi classes must be nonnegative, correct ", a))
        end
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    psi_deg = sum(a)==0 ? 0 : 1
    return Cycle(sum(a), psi_deg)
end

function Jet(n::Int64, deg::Int64, n_marks::Int64, c::Int64, s::Int64, p::Int64, q::Int64)::Cycle
    
    try 
        1 > n_marks && error("Jet is defined when there is at least one mark")
        0 > p       && error(string("exponents of Jet must be nonnegative, correct ", p))
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    psi_deg = p==0 ? 0 : 1
    return Cycle(p+1, psi_deg)
end