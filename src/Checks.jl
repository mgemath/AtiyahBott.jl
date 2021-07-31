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
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);
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
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);
julia> is_zero_cycle(4,1,0,P)
true
```
"""
function is_zero_cycle(n::Int64, deg::Int64, n_marks::Int64, P)::Bool
    
    if !isa(P, Array)  #we want that P is an array, possibly with only one element
        P = [P]
    end

    for res in 1:length(P)
        local P_cycle = P[res](n, deg, n_marks, 0, 0)

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

"""
    check_Data()

List of all files containing the colorations in the folder Data.
"""
function check_Data()::Nothing

    data_dir::String = dirname(current_dir)*"/Data/" #path of the folder Data
    
    if !isdir(data_dir)
        println("Folder ", data_dir, " not found.")
        return
    end
    
    Dimension_dirs = [x for x in readdir(data_dir) if startswith(x,"Dimension")]
    
    if length(Dimension_dirs) == 0
        println("""Folder "Data" is empty.""")
        return
    end
    
    println("""Folder "Data" found.""")
    
    for current_dir in Dimension_dirs
        println(current_dir," contains:")
        files = [x for x in readdir(data_dir*current_dir)]
        for v in 2:14
            num = count(x->length(x)==v+3 && endswith(x,".gz"), files)
            if num == 0
                println( "No colored graphs with ", v, " vertices.")
                continue
            elseif num < number_trees[v-1]
                println("Some colored graph with ",v," vertices is missing.")
            else
                println("All colored graph with ",v," vertices.")
            end
            
        end
    end
    
    return nothing
end

"""
    fill_Data(n, d)

Download from internet all colorations used for computations in the moduli space with dimension `n` and degree `d`.
Return `true` if there is no need to download any file or if all downloads go well. Otherwise, return `false`.
"""
function fill_Data(n::Int64, d::Int64)::Bool

    if n < 1
        printstyled("ERROR: ", bold=true, color=:red)
        println("n must be positive, correct ", n)
        return false
    end
    if d > 13 || d < 1
        printstyled("ERROR: ", bold=true, color=:red)
        println("d must be between 1 and 13, correct ", d)
        return false
    end

    local list_miss::Vector{String} = String[]
    local final_state::Bool = true
    
    local Dim_dir::String = dirname(current_dir)*"/Data/Dimension$n" #path of the folder containing the colorations
    mkpath(Dim_dir) #create the folder

    list_g::IOStream = open(current_dir*"/list_trees.txt") 
    #open the file containing the list of Prufer sequences of graphs
    for v in 2:(d+1) #run the computation among all graphs with fixed number of vertices
        for _ in 1:number_trees[v - 1]  #run the computation for a fixed graph
            str = readline(list_g) #read a new line, we expect a Prufer seq plus the number of automorphisms
            name_file = string(split(str, ',')[1],"0.gz")
            if !(name_file in readdir(Dim_dir))
                push!(list_miss, name_file)
            end
        end
    end

    close(list_g)

    if !isempty(list_miss)
        prog = ProgressUnknown("Downloading colorations...", spinner=true, color=:white)
        for name_file in list_miss
            next!(prog)
            url = "https://raw.githubusercontent.com/mgemath/Colorations/main/Dimension$n/$name_file"
            dest = Dim_dir*"/$name_file"
            try
                Downloads.download(url, dest)
            catch e
                finish!(prog, desc = "Download failed           ", spinner = '✗')
                printstyled(stderr,"ERROR: ", bold=true, color=:red)
                printstyled(stderr,sprint(showerror,e), color=:light_red)
                println(stderr)
                final_state = false
                break  #end for n_g
            end
            final_state || break
        end

        if final_state
            finish!(prog, desc = "All missing colorations have been downloaded.")
        end
    end

    return final_state
end


"""
    free_Data()

Delete the folder Data.
"""
function free_Data()::Nothing
    rm(dirname(current_dir)*"/Data/", force=true, recursive=true)
    return nothing
end
