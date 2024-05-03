export 
    Incidency,
    Hypersurface,
    Contact,
    O1_i,
    O1,
    R1,
    Psi,
    Jet,
    Euler_inv

using Graphs: SimpleEdge

import Base: *, //, /, ^, +, -, inv, one, zero

struct EquivariantClass
    rule::Expr
    func::Function
end 

function Base.show( io::IO, ec::EquivariantClass )::Nothing  
        #print(io, "EquivariantClass[ ... ]" )
        return nothing
end


#########################################
### Operations of equivariant classes ###
#########################################

### Products
function *( ec1::EquivariantClass, ec2::EquivariantClass )::EquivariantClass
    
    rule = quote
        $(ec1.rule)*$(ec2.rule)
    end 

    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end

function *( ec1::EquivariantClass, n::Number )::EquivariantClass
    
    rule = quote
        $(ec1.rule)*$(n)
    end 

    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end

function *( n::Number, ec1::EquivariantClass )::EquivariantClass

    return ec1*n
end

### Division
function //( ec1::EquivariantClass, ec2::EquivariantClass )::EquivariantClass
    
    rule = quote
        $(ec1.rule)//$(ec2.rule)
    end 

    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end

function //( ec1::EquivariantClass, n::Number )::EquivariantClass
    
    rule = quote
        $(ec1.rule)//$(n)
    end 

    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end

### Sums
function +( ec1::EquivariantClass, ec2::EquivariantClass )::EquivariantClass
    
    rule = quote
        $(ec1.rule)+$(ec2.rule)
    end 

    return EquivariantClass( rule, eval(:(( g, c, w, s, m ) -> $rule )))
end

function +( ec1::EquivariantClass, n::Number )::EquivariantClass #it makes sene only if n==0
    
    rule = quote
        $(ec1.rule)+$(n)
    end 

    return EquivariantClass( rule, eval(:(( g, c, w, s, m ) -> $rule )))
end
function +( n::Number, ec1::EquivariantClass )::EquivariantClass #it makes sene only if n==0

    return ec1+n
end
### Minus
function -( ec1::EquivariantClass, ec2::EquivariantClass )::EquivariantClass
    
    rule = quote
        $(ec1.rule)-$(ec2.rule)
    end 

    return EquivariantClass( rule, eval(:(( g, c, w, s, m ) -> $rule )))
end
function -( ec1::EquivariantClass, n::Number )::EquivariantClass #it makes sene only if n==0
    
    rule = quote
        $(ec1.rule)-$(n)
    end 

    return EquivariantClass( rule, eval(:(( g, c, w, s, m ) -> $rule )))
end
function -( n::Number, ec1::EquivariantClass )::EquivariantClass #it makes sene only if n==0

    rule = quote
        $(n)-$(ec1.rule)
    end 

    return EquivariantClass( rule, eval(:(( g, c, w, s, m ) -> $rule )))
end
### Exponent
function ^( ec1::EquivariantClass, n::Number )::EquivariantClass
    
    rule = quote
        $(ec1.rule)^$(n)
    end 

    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end

### Constants
function one(ec1::EquivariantClass)::EquivariantClass
    
    rule = quote
        1
    end 

    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end
function inv(ec1::EquivariantClass)::EquivariantClass
    
    rule = quote
        $(ec1.rule)^(-1)
    end 

    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end
function zero(ec1::EquivariantClass)::EquivariantClass
    
    rule = quote
        0
    end 

    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end
#######################################
### List of all equivariant classes ###
#######################################
"""
    Incidency(r)

Equivariant class of the cycle parameterizing curves meeting a linear subspace of codimension `r`.
# Arguments
- `r::Int64`: the codimension of the subvariety. Alternatively, it can be an array of integers, meaning the multiplication of the equivariant class defined by each element of the array.

# Example
The following Gromov-Witten invariants
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{3},1)}\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{3})^{2} &= 1 \\\\
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{3},1)}\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{2})^{2}\\cdot \\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{3}) &= 1 \\\\
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{3},3)}\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{2})^{2}\\cdot \\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(3))) &= 756 \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using AtiyahBott)
julia> P = Incidency(3)^2;

julia> AtiyahBottFormula(3, 1, 0, P, show_bar=false);
Result: 1

julia> P = Incidency([2,2,3]);

julia> AtiyahBottFormula(3, 1, 0, P, show_bar=false);
Result: 1

julia> P = Incidency([2,2])*Hypersurface(3);

julia> AtiyahBottFormula(3, 3, 0, P, show_bar=false);
Result: 756
```
!!! warning "Attention!"

    The program will stop if `r` is not positive.

"""
function Incidency( r )::EquivariantClass
    
    rule = :(Incidency( g, c, w, s, $r ))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end
function Incidency(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, r::Int64)::fmpq
    
    local p1::fmpq = fmpq(0); #the final result
    local temp1::fmpq = fmpq(1)
    local temp2::fmpq = fmpq(1)
    r -= 1                          
    
    # col = Dict(vertices(g).=> col) #assign colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges

    for e in edges(g)
        for t in (0:r)
            eq!(temp1, scalars[col[src(e)]])
            eq!(temp2, scalars[col[dst(e)]])
            pow_eq!(temp1, t)
            pow_eq!(temp2, r-t)
            mul_eq!(temp1, temp2)
            mul_eq!(temp1, d[e])
            add_eq!(p1, temp1)
            # p1 += d[e]*(scalars[col[src(e)]]^(t))*(scalars[col[dst(e)]]^(r-t))
        end
    end

    return p1

end

function Incidency(g::SimpleGraph{Int64},col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, r::Vector{Int64})::fmpq
    
    local p1::fmpq = fmpq(1)
    local temp1::fmpq = fmpq(1)

    for j in unique(r)
        eq!(temp1, Incidency(g, col, weights, scalars, j))
        pow_eq!(temp1, count(x->x==j, r))
        mul_eq!(p1, temp1)
        # p1 *= Hypersurface(g, col, weights, scalars, j)^count(x->x==j, b)
    end
    return p1

end


"""
    Hypersurface(b)

Equivariant class of the Euler class of the bundle equal to the direct image under the forgetful map of ``\\mathrm{ev}^*\\mathcal{O}_{\\mathbb{P}^n}(b)``. It parameterizes curves contained in a hypersurface of degree `b`.
# Arguments
- `b::Int64`: the degrees of the hypersurface. Alternatively, it can be an array of integers, meaning the multiplication of the equivariant class defined by each element of the array.

# Example
The following Gromov-Witten invariants of Calabi-Yau threefolds
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{4},1)}\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{4}}(5))) &= 2875 \\\\
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{5},2)}\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{5}}(3))^{\\oplus 2}) &= \\frac{423549}{8} \\\\
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{5},3)}\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{5}}(4)))\\cdot\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{5}}(2))) &= \\frac{422690816}{27} \\\\
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{7},4)}\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{7}}(2)))^4 &= 25705160 \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using AtiyahBott)
julia> P = Hypersurface(5);

julia> AtiyahBottFormula(4, 1, 0, P, show_bar=false);
Result: 2875

julia> P = Hypersurface([3,3]);

julia> AtiyahBottFormula(5, 2, 0, P, show_bar=false);
Result: 423549//8

julia> P = Hypersurface(4)*Hypersurface(2);

julia> AtiyahBottFormula(5, 3, 0, P, show_bar=false);
Result: 422690816//27

julia> P = Hypersurface(2)^4;

julia> AtiyahBottFormula(7, 4, 0, P, show_bar=false);
Result: 25705160
```
!!! warning "Attention!"

    The program will stop if `b` is not positive.

"""
function Hypersurface( b )::EquivariantClass
    
    rule = :(Hypersurface( g, c, w, s, $b ))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end

function Hypersurface(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, b::Int64)::fmpq

    local p1::fmpq = fmpq(1)
    # local q1::fmpq = fmpq(1)
    local temp1::fmpq = fmpq(1)
    local temp2::fmpq = fmpq(1)
    
    # col = Dict(vertices(g).=> col) #assign colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    for e in edges(g)
        for alph in 0:(b*d[e])
            eq!(temp1, scalars[col[src(e)]])
            eq!(temp2, scalars[col[dst(e)]])
            mul_eq!(temp1, alph)
            mul_eq!(temp2, b*d[e]-alph)
            add_eq!(temp1, temp2)
            div_eq!(temp1, d[e])
            mul_eq!(p1, temp1)
            # p1 *= (alph*scalars[col[src(e)]]+(b*d[e]-alph)*scalars[col[dst(e)]])//d[e]
        end
    end
    
    for v in vertices(g)
        eq!(temp1, scalars[col[v]])
        mul_eq!(temp1, b)
        pow_eq!(temp1, 1-length(all_neighbors(g, v)))
        mul_eq!(p1, temp1)
        #q1 *= (b*scalars[col[v]])^(1-length(all_neighbors(g, v)))   
    end

    # mul_eq!(p1,q1)
    
    return p1

end

function Hypersurface(g::SimpleGraph{Int64},col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, b::Vector{Int64})::fmpq

    local p1::fmpq = fmpq(1)
    local temp1::fmpq = fmpq(1)

    for j in unique(b)
        eq!(temp1, Hypersurface(g, col, weights, scalars, j))
        pow_eq!(temp1, count(x->x==j, b))
        mul_eq!(p1, temp1)
        # p1 *= Hypersurface(g, col, weights, scalars, j)^count(x->x==j, b)
    end
    return p1
end


"""
    Contact()

Equivariant class of the Euler class of the bundle equal to the direct image under the forgetful map of ``\\mathrm{ev}^*\\mathcal{O}_{\\mathbb{P}^n}(2)`` tensor the dualizing sheaf of the forgetful map. It parameterizes contact curves in an odd dimensional projective space.

# Example
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^{3},1)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{2}\\cdot\\mathrm{ev}_{2}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{3}\\cdot\\mathrm{c_{top}}(\\delta_{*}(\\omega_{\\delta}\\otimes\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(2))) &= 1 \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using AtiyahBott)
julia> P = O1_i(1)^2*O1_i(2)^3*Contact();

julia> AtiyahBottFormula(3, 1, 2, P, show_bar=false);
Result: 1
```
"""
function Contact()::EquivariantClass
    
    rule = :(Contact( g, c, w, s))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end
function Contact(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}})::fmpq

    local p1::fmpq = fmpq(1)
    # local q1::fmpq = fmpq(1)
    local temp1::fmpq = fmpq(1)
    local temp2::fmpq = fmpq(1)
    
    # col = Dict(vertices(g).=> col) #assign colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    for e in edges(g)
        for alph in (1:2*d[e]-1)
            eq!(temp1, scalars[col[src(e)]])
            eq!(temp2, scalars[col[dst(e)]])
            mul_eq!(temp1, alph)
            mul_eq!(temp2, 2*d[e]-alph)
            add_eq!(temp1, temp2)
            div_eq!(temp1, d[e])
            mul_eq!(p1, temp1)
            # p1 *= (alph*scalars[col[src(e)]]+(2*d[e]-alph)*scalars[col[dst(e)]])//d[e]
        end
    end
    
    for v in vertices(g)
        eq!(temp1, scalars[col[v]])
        mul_eq!(temp1, 2)
        pow_eq!(temp1, length(all_neighbors(g, v))-1)
        mul_eq!(p1, temp1)
        # q1 *= (2*scalars[col[v]])^(length(all_neighbors(g, v))-1)
    end

    # mul_eq!(p1, q1)
    return p1

end

"""
    O1_i(j)

Equivariant class of the pull-back of ``\\mathcal{O}_{\\mathbb{P}^n}(1)`` with respect to the j-th evaluation map.
# Arguments
- `j::Int64`: the evaluation map.

# Example
The following Gromov-Witten invariants
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^{1},1)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{1}}(1)\\cdot\\mathrm{ev}_{2}^{*}\\mathcal{O}_{\\mathbb{P}^{1}}(1) &= 1 \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{3},1)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2}\\cdot\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(2))) &= 4
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using AtiyahBott)
julia> P = O1_i(1)*O1_i(2);

julia> AtiyahBottFormula(1, 1, 2, P, show_bar=false);
Result: 1

julia> P = O1_i(1)^2*Hypersurface(2);

julia> AtiyahBottFormula(3, 1, 1, P, show_bar=false);
Result: 4
```
!!! warning "Attention!"

    The program will stop if `j` is not between 1 and the number of marks.

"""
function O1_i( j::Int64)::EquivariantClass
    rule = :(O1_i(g, c, w, s, m, $j ))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end

function O1_i(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, mark::Marks_type, j::Int64)::fmpq
#function O1_i(col::Vector{UInt8}, scalars::Tuple{Vararg{fmpq}}, mark::Marks_type, j::Int64)::fmpq    
    return scalars[col[mark[j]]]
end


"""
    O1()

Equivariant class of the pull-back of ``\\mathcal{O}_{\\mathbb{P}^n}(1)`` with respect to the product of all evaluation maps.

This function is equivalent to the product of the function `O1_i(j)` where `j` runs from 1 to the number of marks.

# Example
The following Gromov-Witten invariants
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,8}(\\mathbb{P}^{2},3)}\\prod_{i=1}^{8}\\mathrm{ev}_{i}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^2 &= 12 \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{3},2)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^2\\cdot\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(3))) &= 81 \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using AtiyahBott)
julia> P = O1()^2;

julia> AtiyahBottFormula(2, 3, 8, P, show_bar=false);
Result: 12

julia> P = O1()^2*Hypersurface(3);

julia> AtiyahBottFormula(3, 2, 1, P, show_bar=false);
Result: 81
```

In order to remove `O1_i(j)` for some `j`, it is enough to divide by that function.

# Example
```julia-repl
julia> P = O1()//O1_i(1);
```
Here `P` is the product of all `O1_i(j)` where `j` runs from 2 to `m`.
"""
function O1()::EquivariantClass
    rule = :(O1( g, c, w, s, m ))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end

function O1(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, mark::Marks_type)::fmpq
    
    local p1 = fmpq(1)
    # col = Dict(vertices(g).=> col)
    
    for t in 1:length(mark)
        mul_eq!(p1, scalars[col[mark[t]]])
        # p1 *= scalars[col[mark[t]]]
    end
    
    return p1
end

"""
    R1(k)

The equivariant class of the first derived functor of the pull-back of ``\\mathcal{O}_{\\mathbb{P}^n}(-k)``.
# Arguments
- `k::Int64`: a positive integer.

# Example
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{1},d)}\\mathrm{c_{top}}(R^{1}\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(-1)))^2 &= \\frac{1}{d^3} \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using AtiyahBott)
julia> d = 1; #for other values of d, change this line

julia> P = R1(1)^2;

julia> AtiyahBottFormula(1, d, 0, P, show_bar=false);
Result: 1
```
!!! warning "Attention!"

    The program will stop if `k` is not positive.

"""
function R1( k )::EquivariantClass
    
    rule = :(R1( g, c, w, s, $k ))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end
function R1(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, k::Int64)::fmpq
    
    local p1 = fmpq(1)
    # local q1 = fmpq(1)
    local temp1::fmpq = fmpq(1)
    local temp2::fmpq = fmpq(1)
    
    # col = Dict(vertices(g).=> col) #assign colors to vertices
    d = Dict(edges(g).=> weights) #assign weights to edges
    
    for e in edges(g)
        for alph in 1:(k*d[e]-1)
            eq!(temp1, scalars[col[src(e)]])
            eq!(temp2, scalars[col[dst(e)]])
            mul_eq!(temp1, alph)
            mul_eq!(temp2, k*d[e]-alph)
            add_eq!(temp1, temp2)
            div_eq!(temp1, d[e])
            neg!(temp1)
            mul_eq!(p1, temp1)
            # p1 *= (-1)*(alph*scalars[col[src(e)]]+(k*d[e]-alph)*scalars[col[dst(e)]])//d[e]
        end
    end
    
    for v in vertices(g)
        eq!(temp1, scalars[col[v]])
        mul_eq!(temp1, -k)
        pow_eq!(temp1, length(all_neighbors(g, v))-1)
        mul_eq!(p1, temp1)
        # q1 *= (-k*scalars[col[v]])^(length(all_neighbors(g, v))-1)   
    end

    # mul_eq!(p1, q1)
    return p1

end

"""
    Psi(a)

Equivariant class of the cycle of ``\\psi``-classes.
# Arguments
- `a::Vector{Int64}`: the vector of the exponents of the ``\\psi`` classes. It is ordered, meaning that the first element is the exponent of ``\\psi_1``, the second is the exponent of ``\\psi_2``, and so on.

!!! note

    The size of `a` must be at most `m`. If it is smaller, missing exponents will be considered as zeros.
    If `a` is a number, it will be considered as the exponent of ``\\psi_1``.

!!! warning "Attention!"

    The program will stop if we have one of the following conditions:

    * the size of `a` is bigger than `m`,
    * `a` contains a negative number.

# Example
The following Gromov-Witten invariants
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^{6},2)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{6}}(1)^{5}\\cdot\\mathrm{ev}_{2}^{*}\\mathcal{O}_{\\mathbb{P}^{6}}(1)^{2}\\cdot\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{6}}(5)))\\cdot\\psi_{1}\\psi_{2}^{0} &= 495000 \\\\
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^{10},2)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{10}}(1)^{8}\\cdot\\mathrm{ev}_{2}^{*}\\mathcal{O}_{\\mathbb{P}^{10}}(1)^{6}\\cdot\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{10}}(7)))\\cdot\\psi_{1}^{2} &= 71804533752 \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{2},2)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2}\\cdot\\psi_{1}^{4} &= \\frac{1}{8} \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{2},2)}\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2})^{4}\\cdot\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)\\cdot(\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)+\\psi_{1}) &= 2 \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{3},2)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)\\cdot(\\psi_{1}^{7}\\cdot\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)+\\psi_{1}^{6}\\cdot\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2}) &= -\\frac{5}{16} \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using AtiyahBott)
julia> P = O1_i(1)^5*O1_i(2)^2*Hypersurface(5)*Psi([1,0]);

julia> AtiyahBottFormula(6, 2, 2, P, show_bar=false);
Result: 495000

julia> P = O1_i(1)^8*O1_i(2)^6*Hypersurface(7)*Psi(2);

julia> AtiyahBottFormula(10, 2, 2, P, show_bar=false);
Result: 71804533752

julia> P = O1()^2*Psi(4);

julia> AtiyahBottFormula(2, 2, 1, P, show_bar=false);
Result: 1//8

julia> P = Incidency(2)^4*O1_i(1)*(O1_i(1) + Psi(1));

julia> AtiyahBottFormula(2, 2, 1, P, show_bar=false); #number of plane conics through four points and tangent to a line
Result: 2

julia> P = O1()*(Psi(7)*O1()+Psi(6)*O1()^2);

julia> AtiyahBottFormula(3, 2, 1, P, show_bar=false);
Result: -5//16
```
!!! warning "Psi is singleton!"

    `Psi` cannot be multiplied by itself.
    ```julia-repl
    julia> P = O1()^2*Psi(1)^4;                  #this is **wrong**
    julia> AtiyahBottFormula(2,2,1,P);
    Warning: more instances of Psi has been found. Type:
    julia> ?Psi
    for support.
    julia> P = O1()^2*Psi(3)*Psi(1);             #this is **wrong**
    julia> AtiyahBottFormula(2,2,1,P);
    Warning: more instances of Psi has been found. Type:
    julia> ?Psi
    for support.
    julia> P = O1()^2*Psi(4);
    julia> AtiyahBottFormula(2,2,1,P);
    Result: 1//8
    ```
"""
function Psi( a )::EquivariantClass
    
    rule = :(Psi( g, c, w, s, m, $a ))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end
function Psi(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, mark::Marks_type, a::Int64)::fmpq
    
    return Psi(g, col, weights, scalars, mark, [a])
end
function Psi(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, mark::Marks_type, a::Vector{Int64})::fmpq
        
    findfirst(x -> x>0, a) === nothing && return fmpq(1) #if all of them are zero or a is empty
    
    local q1::fmpq = fmpq(1)
    local temp1::fmpq = fmpq(1)
    local Sum_ai::Int64
    local n::Int64
    local M::Int64
    local d = Dict(edges(g).=> weights) #assign weights to edges
    local inv_marks::Dict{Int64,Vector{Int64}} = invert_marks(mark, nv(g))
    
    for v in vertices(g)
        
        a_v = Int64[]
        for i in inv_marks[v]
            (i > length(a) || a[i] == 0) && continue
            push!(a_v, a[i])
        end

        Sum_ai = sum(a_v)
        Sum_ai == 0 && continue #if S contains only zeros, or it is empty, continue

        n = length(all_neighbors(g, v)) + length(inv_marks[v])
        
        n > 2 && Sum_ai > n - 3 && return fmpq(0)
        
        #If no previous condition holds, then n>1
        if n == 2 #necessary |S_v| == 1
            M = (-1)^a_v[1]
        else # n>2 and Sum_ai <= n - 3
            M = multinomial(n - 3 - Sum_ai, a_v...)
        end
        

        local s1 = fmpq(0)
        
        for w in all_neighbors(g, v)
            # e = SimpleEdge(v,w)
            # d_e = haskey(d,e) ? d[e] : d[reverse(e)]
            eq!(temp1, scalars[col[w]])
            neg!(temp1)
            add_eq!(temp1, scalars[col[v]])
            div_eq!(temp1, d[SimpleEdge(min(v,w),max(v,w))])
            inv!(temp1)
            add_eq!(s1, temp1)
            # s1 += d_e//(scalars[col[v]]-scalars[col[w]])
        end
        pow_eq!(s1, -Sum_ai)
        mul_eq!(q1, s1)
        mul_eq!(q1, M)
        # s1 ^= -Sum_ai
        # q1 *= M*s1

    end
        
    return q1
end

"""
    Jet(p, q)

Equivariant class of the jet bundle ``J^p`` of the pull back of ``\\mathcal{O}_{\\mathbb{P}^n}(q)`` with respect to the first ``\\psi``-class.
# Arguments
- `p::Int64`: the exponent of the Jet bundle. In particular, it is a bundle of rank ``p+1``.
- `q::Int64`: the degree of the line bundle that is pulled back.


!!! note

    In order to define this bundle, the number of marks must be at least 1.
    You cannot multiply this bundle by the class `Psi(a)`.


# Example
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{2},2)}\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2})^{4}\\cdot\\mathrm{c_{top}}(J^{1}(\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1))) &= 2 \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{2},2)}\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2})^{4}\\cdot(\\mathrm{c_{top}}(J^{1}(\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)))+\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2}) &= 3 \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{3},d)}\\frac{\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{2}}{k}\\cdot\\mathrm{c_{top}}(J^{4d-2}(\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(k))) &= \\frac{(4d-2)!}{(d!)^{4}} \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using AtiyahBott)
julia> P = Incidency(2)^4*Jet(1,1);

julia> AtiyahBottFormula(2, 2, 1, P, show_bar=false);
Result: 2

julia> P = Incidency(2)^4*(Jet(1,1)+O1()^2);

julia> AtiyahBottFormula(2, 2, 1, P, show_bar=false);
Result: 3

julia> d=1;k=1; #for other values of d, change this line

julia> P = (O1()^2)//k*Jet(4*d-2,k);

julia> AtiyahBottFormula(3, d, 1, P, show_bar=false);   #The value of this integral does not depend on k, only on d
Result: 2
```
"""
function Jet( p, q )::EquivariantClass
    
    rule = :(Jet( g, c, w, s, m, $p, $q ))
    return EquivariantClass( rule, eval( :(( g, c, w, s, m ) -> $rule )))
end
function Jet(g::SimpleGraph{Int64}, col::Tuple{Vararg{Int64}}, weights::Vector{Int64}, scalars::Tuple{Vararg{fmpq}}, mark::Marks_type, p::Int64, q::Int64)::fmpq
    
    local s1::fmpq = fmpq(0)
    local temp1::fmpq = fmpq(1)

    for h in 0:p
        eq!(temp1, O1_i(g, col, weights, scalars, mark, 1))
        mul_eq!(temp1, q)
        pow_eq!(temp1, p+1-h)
        mul_eq!(temp1, stirlings1(p+1, p+1-h))
        mul_eq!(temp1, Psi(g, col, weights, scalars, mark, [h]))
        add_eq!(s1, temp1)
        # s1 += stirlings1(p+1, p+1-h)*(q*O1_i(g, col, weights, scalars, mark, 1))^(p+1-h)*Psi(g, col, weights, scalars, mark, [h])
    end

    return s1
end




###############################################################
### List of the class useful for computing the rank (Cycle) ###
###############################################################
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