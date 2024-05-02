"""
**AtiyahBott** is a module containing an implementation of the Atiyah-Bott residue formula in the Julia language.
The theory and the algorithm behind the package is described in the paper 
"Effective computations of the Atiyah-Bott formula" by Giosuè Muratore and 
Csaba Schneider (https://doi.org/10.1016/j.jsc.2022.01.005).

A previous version used colorations from here: https://github.com/mgemath/Colorations. The current version generates internally all colorations.
"""
module AtiyahBott

using Combinatorics 
using Graphs 
using ProgressMeter
using Nemo

include("Arithmetic.jl")
include("Marked.jl")
include("GraphFunctions.jl")
include("EquivariantClasses.jl")
include("Checks.jl")
include("Main.jl")
include("NemoFunctions.jl")
include("Trees.jl")
include("Colors.jl")


end # module