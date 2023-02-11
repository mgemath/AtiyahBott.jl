# AtiyahBott.jl
[![Doc](https://img.shields.io/badge/docs-stable-blue.svg)](https://mgemath.github.io/AtiyahBott.jl/)
[![Site](https://img.shields.io/website?up_message=Arxiv&url=https%3A%2F%2Farxiv.org%2Fabs%2F2105.11183)](https://arxiv.org/abs/2105.11183)
[![Site](https://img.shields.io/website?up_message=JSC&url=https%3A%2F%2Fwww.sciencedirect.com%2Fscience%2Farticle%2Fpii%2FS0747717122000050)](https://www.sciencedirect.com/science/article/pii/S0747717122000050)

This package contains an implementation of the Atiyah-Bott residue formula for the moduli space of genus 0 stable maps in the Julia language. The theory behind the package and the algorithm are described in the paper 
"Effective computations of the Atiyah-Bott formula" by Giosuè Muratore and Csaba Schneider (https://doi.org/10.1016/j.jsc.2022.01.005).<br>
Full documentation is available here: https://mgemath.github.io/AtiyahBott.jl/.

## Installation
In order to install this package, type:
```julia
julia> using Pkg
julia> Pkg.add("AtiyahBott")
```
After the installation, simply type:
```julia
julia> using AtiyahBott
```
every time you want to use the program.

To use our code, you should first define the equivariant classes to be calculated as 
```julia
julia> P = ...
```
After the "=", one has to write an expression in the equivariant classes. After P is defined, one has to call the
Atiyah-Bott formula by the command
```julia
julia> AtiyahBottFormula(n,d,m,P);
```
The full list of the currently supported equivariant classes is the following:
```julia
O1_i(j)       (pull back of the line bundle O(1) with respect to the ev_j)
O1()           (product of all O1_i(j))
Incidency(r)    (class of curves meeting a linear subspace)
Hypersurface(b) (class of curves contained in a hypersurface of degree b)
Contact()        (class of contact curves)
R1(k)           (first derived functor of direct image of the pull back of O(-k))
Psi(a)        (cycle of psi-classes)
Jet(p,q)      (Euler class of the jet bundle J^p)
```
Brief descriptions on these functions can be obtained through the standard help functionality of Julia by typing "?" and then the name of the function.
```julia
help?> Psi
```

Note that computations can be faster using multi-threading. Visit https://docs.julialang.org/en/v1/manual/multi-threading/#man-multithreading to learn how to start Julia with multi-threading.
## Examples
In the following we list some geometrically meaning computations.

### Curves in projective spaces

To compute the number of rational plane curves of degree d through 3d−1 general points, one may write:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = O1()^2;
julia> AtiyahBottFormula(2,d,3*d-1,P);
```
Alternatively, one can perform such computation with zero marked points by typing:
```julia
julia> P = Incidency(2)^(3*d-1);
julia> AtiyahBottFormula(2,d,0,P);
```
### Curves in Hypersurfaces

The virtual number of rational degree d curves on a general complete intersection of type (2,3) in the projective space of dimension 5:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = Hypersurface([2,3]);
julia> AtiyahBottFormula(5,d,0,P);
```
The number of rational degree d curves on a cubic surface passing through d-1 points:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = Hypersurface(3)*(Incidency(2)//3)^(d-1);
julia> AtiyahBottFormula(3,d,0,P);
```

### Tangency conditions

The number plane rational degree d curves through 3d-2 points and tangent to a line:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = Incidency(2)^(3*d-1)*Jet(1,1);
julia> AtiyahBottFormula(2,d,1,P);
```

### Hurwitz numbers
The weighted number of genus 0 degree d covers of the projective line, which are étale over a fixed point and with 2d-2 fixed finite simple ramification points, is:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = O1()*Psi(ones(Int,2*d-2));
julia> AtiyahBottFormula(1,d,2*d-2,P);
```
See https://arxiv.org/pdf/math/0101147.pdf.

## Future goals
The following may be future expansions of this program.
 - Support for positive genus curves.
 - Improve parallel acceleration.

If you have other suggestions, please raise an issue on github. 

# Citing

We encourage you to cite our work if you have used our package. See "Cite this repository" on this page.
