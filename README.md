# AtiyahBott.jl
[![Doc](https://img.shields.io/badge/docs-stable-blue.svg)](https://mgemath.github.io/AtiyahBott.jl/)
[![Site](https://img.shields.io/website?up_message=Arxiv&url=https%3A%2F%2Farxiv.org%2Fabs%2F2105.11183)](https://arxiv.org/abs/2105.11183)

This package contains an implementation of the Atiyah-Bott residue formula for the moduli space of genus 0 stable maps in the Julia language. The theory behind the package and the algorithm are described in the paper 
"Effective computations of the Atiyah-Bott formula" by Giosuè Muratore and Csaba Schneider (https://arxiv.org/pdf/2105.11183.pdf).<br>
Full documentation is available here: https://mgemath.github.io/AtiyahBott.jl/.

## Installation
In order to install this package, type:
```julia
julia> using Pkg
julia> Pkg.add("AtiyahBott")
```

The installation will download automatically all colorations up to degree 4 and dimension 4 from here: https://github.com/mgemath/Colorations.

After the installation, simply type:
```julia
julia> using AtiyahBott
```
every time you want to use the program.

To use our code, you should first define the equivariant classes to be calculated as 
```julia
julia> P = (g,c,w,s,m) ->...
```
After the "->", one has to write an expression in the equivariant classes. After P is defined, one has to call the
Atiyah-Bott formula by the command
```julia
julia> AtiyahBottFormula(n,d,m,P);
```
The full list of the currently supported equivariant classes is the following:
```julia
O1_i(g,c,w,s,m,i)       (pull back of the line bundle O(1) with respect to the ev_i)
O1(g,c,w,s,m)           (product of all O1_i(g,c,w,s,m,i))
Incidency(g,c,w,s,r)    (class of curves meeting a linear subspace)
Hypersurface(g,c,w,s,b) (class of curves contained in a hypersurface)
Contact(g,c,w,s)        (class of contact curves)
R1(g,c,w,s,k)           (first derived functor of direct image of the pull back of O(-k))
Psi(g,c,w,s,m,a)        (cycle of psi-classes)
Jet(g,c,w,s,m,p,q)      (Euler class of the jet bundle J^p)
```
Brief descriptions on these functions can be obtained through the standard help functionality of Julia by typing "?" and then the name of the function.
```julia
help?> Psi
```
## Examples
In the following we list some geometrically meaning computations.

### Curves in projective spaces

To compute the number of rational plane curves of degree d through 3d−1 general points, one may write:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2;
julia> AtiyahBottFormula(2,d,3*d-1,P);
```
Alternatively, one can perform such computation with zero marked points by typing:
```julia
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^(3*d-1);
julia> AtiyahBottFormula(2,d,0,P);
```
### Curves in Hypersurfaces

The virtual number of rational degree d curves on a general complete intersection of type (2,3) in the projective space of dimension 5:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,[2,3]);
julia> AtiyahBottFormula(5,d,0,P);
```
The number of rational degree d curves on a cubic surface passing through d-1 points:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,3)*(Incidency(g,c,w,s,2)//3)^(d-1);
julia> AtiyahBottFormula(3,d,0,P);
```

### Tangency conditions

The number plane rational degree d curves through 3d-2 points and tangent to a line:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^(3*d-1)*Jet(g,c,w,s,m,1,1);
julia> AtiyahBottFormula(2,d,1,P);
```

### Hurwitz numbers
The weighted number of genus 0 degree d covers of the projective line, which are étale over a fixed point and with 2d-2 fixed finite simple ramification points, is:
```julia
julia> d = 1; #for other values of d, change this line
julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)*Psi(g,c,w,s,m,ones(Int,2*d-2));
julia> AtiyahBottFormula(1,d,2*d-2,P);
```
See https://arxiv.org/pdf/math/0101147.pdf.

## Future goals
The following may be future expansions of this program.
 - Support for positive genus curves.
 - GPU or parallel acceleration.
 - Alternative use of the type "Rational{BigInt}", which is a bottleneck.

If you have other suggestions, please raise an issue on github. 
