export
    dim_M,
    codim,
    is_zero_cycle,
    check_Data,
    fill_Data,
    free_Data

import Base: *, //, /, ^, +, -, inv, one, zero

struct Cycle #We define the structure Cycle. It keeps track of the codimension and of the fact that it contains psi-classes
    c::Int64 #The codimension. If negative, the class is not well defined
    is_psi::Int64 #this is a counter that measure how many times the class Psi is called in P. If this number is strictly greater than 1, then the program stops.
    valid::Bool #check if it is a valid cycle
end

function error_cycle()::Cycle
    return Cycle(0, 0, false)
end

function Cycle(c::Int64, is_psi::Int64)::Cycle
    return Cycle(c, is_psi, true)
end

#Cycles satisfy the following arithmetic. They are not invertible (except in codimension 0), and cycles of different codimension are not summable.
#So the Cycle(-1, 0) is a meaningless cycle.

*(P1::Cycle, P2::Cycle)::Cycle = Cycle(P1.c + P2.c, P1.is_psi + P2.is_psi, P1.valid && P2.valid)
^(P1::Cycle, n::Int64) ::Cycle = Cycle(n*P1.c, n*P1.is_psi, P1.valid) #if n=0, we got a zero cycle
+(P1::Cycle, P2::Cycle)::Cycle = P1.c == P2.c ? Cycle(P1.c, max(P1.is_psi, P2.is_psi), P1.valid && P2.valid) : error_cycle()
-(P1::Cycle, P2::Cycle)::Cycle = +(P1, P2)
*(P1::Cycle, ::Number) ::Cycle = P1
*(::Number, P1::Cycle) ::Cycle = P1
//(P1::Cycle, ::Number)::Cycle = P1
//(::Number, P1::Cycle)::Cycle = Cycle(-P1.c, -P1.is_psi, P1.valid)
/(P1::Cycle, ::Number) ::Cycle = P1
/(::Number, P1::Cycle) ::Cycle = Cycle(-P1.c, -P1.is_psi, P1.valid)
+(P1::Cycle, n::Number)::Cycle = n == 0 ? P1 : error_cycle()
-(P1::Cycle, n::Number)::Cycle = +(P1, n)
+(n::Number, P1::Cycle)::Cycle = +(P1, n)
-(n::Number, P1::Cycle)::Cycle = +(P1, n)
# +(P1::Cycle)           ::Cycle = P1
# -(P1::Cycle)           ::Cycle = P1
//(P1::Cycle,P2::Cycle)::Cycle = Cycle(P1.c-P2.c, P1.is_psi-P2.is_psi, P1.valid && P2.valid)
/(P1::Cycle, P2::Cycle)::Cycle = //(P1, P2)
one(::Cycle)           ::Cycle = Cycle(0, 0, true)
inv(P1::Cycle)         ::Cycle = Cycle(-P1.c, -P1.is_psi, P1.valid)
zero(::Cycle)          ::Int64 = 0

"""
    dim_M(n, d, m)

The dimension of the moduli space of stable rational map to the projective space of dimension `n`, of degree `d` with `m` marks.
# Arguments
- `n::Int64`: the dimension of the projective space.
- `d::Int64`: the degree of the stable maps.
- `m::Int64`: the number of marks.

# Example
```julia-repl
julia> dim_M(2,2,5)
10
```
"""
function dim_M(n::Int64, deg::Int64, n_marks::Int64)::Int64
    
    return n + (n + 1)*deg + n_marks - 3
end


"""
    codim(n, d, m, P)

The codimension of the equivariant class `P`.
# Arguments
- `n::Int64`: the dimension of the projective space.
- `d::Int64`: the degree of the stable maps.
- `m::Int64`: the number of marks.
- `P`: the equivariant class.

# Example
```julia-repl
julia> P = Hypersurface(g,c,w,s,5);
julia> codim(4,1,0,P)
6
```
"""
function codim(n::Int64, deg::Int64, n_marks::Int64, P)::Int64
    
    return P(n, deg, n_marks, 0, 0).c
end

"""
    is_zero_cycle(n, d, m, P)

Return `true` if the equivariant class `P` is a 0-cycle in the moduli space, `false` otherwise.
# Arguments
- `n::Int64`: the dimension of the projective space.
- `deg::Int64`: the degree of the stable maps.
- `n_marks::Int64`: the number of marks.
- `P`: the equivariant class.

# Example
```julia-repl
julia> P = Hypersurface(g,c,w,s,5);
julia> is_zero_cycle(4,1,0,P)
true
```
"""
function is_zero_cycle(n::Int64, deg::Int64, n_marks::Int64, P_input)::Bool

    local n_results::Int64 = 0

    if isa(P_input, Array)  #we want that P is an array, possibly with only one element
        n_results = length(P_input)
    else
        n_results = 1
        #P = [P]
    end
    local P::Vector{Function} = Vector(undef, n_results)
    
    if isa(P_input, Array)  #we want that P is an array, possibly with only one element
        if typeof(P_input[1]) == EquivariantClass
            for res in eachindex(P)
                P[res] = P_input[res].func
            end
        else
            P = P_input
        end
    else
        if typeof(P_input) == EquivariantClass
            P[1] = P_input.func
        else
            P[1] = P_input
        end
    end

    for res in eachindex(P)
        local P_cycle = Base.invokelatest( P[res], n, deg, n_marks, 0, 0)

        if !P_cycle.valid
            printstyled("Warning: ", bold=true, color=:light_yellow)
            println("the class is not valid")
            return false
        end

        if P_cycle.is_psi > 1
            printstyled("Warning: ", bold=true, color=:light_yellow)
            println("more instances of Psi has been found. Type:")
            printstyled("julia> ", bold=true, color=:light_green)
            println("?Psi")
            println("for support.")
            return false
        end

        if P_cycle.c != dim_M(n, deg, n_marks)
            printstyled("Warning: ", bold=true, color=:light_yellow)
            length(P)==1 ? println("the class is not a zero cycle") : println("some classes are not zero cycles")
            return false
        end
    end
    return true
end
