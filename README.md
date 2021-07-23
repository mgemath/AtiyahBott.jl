# AtiyahBott.jl
This package contains an implementation of the Atiyah-Bott residue formula for the moduli space of genus 0 stable maps in the Julia language. The theory behind the package is described in the paper 
"Effective computations of the Atiyah-Bott formula" by Giosuè Muratore e Csaba Schneider (https://arxiv.org/pdf/2105.11183.pdf).

In order to install this package, type:
julia> using Pkg<br>
julia> Pkg.add(url="https://github.com/mgemath/AtiyahBott.jl")

After that, simply type

julia> using AtiyahBott;

To use our code, you should first define the equivariant classes to be calculated as "P = (g,c,w,s,m) ->...".<br>
After the "->", one has to write an expression in the equivariant classes. After P is defined, one has to call the
Atiyah-Bott formula by the command AtiyahBottFormula(n,d,m,P). 

The full list of the currently supported equivariant classes is the following:

<pre class="tab">
O1_i(g,c,w,s,m,i)       (pull back of the line bundle O(1) of the projective space)<br>
O1(g,c,w,s,m)           (product of all O1_i)<br>
Incidency(g,c,w,s,r)    (class of curves meeting a linear subspace)<br>
Hypersurface(g,c,w,s,b) (class of curves contained in a hypersurface)<br>
Contact(g,c,w,s)        (class of contact curves)<br>
R1(g,c,w,s,k)           (first derived functor of direct image of the pull back of O(-k))<br>
Psi(g,c,w,s,m,a)        (cycle of psi-classes) <br> 
Jet(g,c,w,s,m,p,q)      (Euler class of the jet bundle J^p)<br>
</pre>
Brief descriptions on these functions can be obtained through the standard help functionality of Julia by typing "?" and then the name of the function.

For example, to compute the number of rational plane curves of degree d through 3d−1 general points, one may write
<pre class="tab">
julia> d = 1 #for other values of d, change this line<br>
julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2<br>
julia> AtiyahBottFormula(2,d,3*d-1,P);<br>
</pre>
Alternatively, one can perform such computation with zero marked points by typing
<pre class="tab">
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^(3*d-1)<br>
julia> AtiyahBottFormula(2,d,0,P);<br>
</pre>
The virtual number of rational degree d curves on a general complete intersection of type (2,3) in the projective space of dimension 5:
<pre class="tab">
julia> d = 1 #for other values of d, change this line<br>
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,[2,3])<br>
julia> AtiyahBottFormula(5,d,0,P);<br>
</pre>
The number of rational degree d curves on a cubic surface passing through d-1 points:
<pre class="tab">
julia> d = 1 #for other values of d, change this line<br>
julia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,3)*(Incidency(g,c,w,s,2)//3)^(d-1)<br>
julia> AtiyahBottFormula(3,d,0,P);<br>
</pre>
The number plane rational degree d curves through 3d-2 points and tangent to a line:
<pre class="tab">
julia> d = 1 #for other values of d, change this line<br>
julia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^(3*d-1)*Jet(g,c,w,s,m,1,1);<br>
julia> AtiyahBottFormula(2,d,1,P);<br>
</pre>

In order to execute these computations, the system requires the list of tree graphs up to d+1 vertices. These graphs are available in the file list_trees.txt. 
We also need the possible colorings of these graphs with n+1 colors. These can either be computed by the package or they can be read from a file stored 
in the Data directory (this is the fastest option). The user can use function check_Data() to see the list of graphs with pre-computed colorings. 
