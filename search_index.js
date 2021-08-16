var documenterSearchIndex = {"docs":
[{"location":"#AtiyahBott.jl","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"","category":"section"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"","category":"page"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"AtiyahBott","category":"page"},{"location":"#AtiyahBott","page":"AtiyahBott.jl","title":"AtiyahBott","text":"AtiyahBott is a module containing an implementation of the Atiyah-Bott residue formula in the Julia language. The theory and the algorithm behind the package is described in the paper  \"Effective computations of the Atiyah-Bott formula\" by Giosuè Muratore and  Csaba Schneider (https://arxiv.org/pdf/2105.11183.pdf).\n\nThe colorations (useful to speed up the code) are up to projective spaces of dimension 29. The full list is here: https://github.com/mgemath/Colorations/.\n\n\n\n\n\n","category":"module"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"In order to install the package, type:","category":"page"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"julia> using Pkg\r\njulia> Pkg.add(\"AtiyahBott\")","category":"page"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"To load the package, type:","category":"page"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"julia> using AtiyahBott","category":"page"},{"location":"#The-function-AtiyahBottFormula","page":"AtiyahBott.jl","title":"The function AtiyahBottFormula","text":"","category":"section"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"This is the main function of the package.","category":"page"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"AtiyahBottFormula(n::Int64, deg::Int64, n_marks::Int64, P, do_check::Bool = true, show_bar::Bool = true, down_col::Bool = true)","category":"page"},{"location":"#AtiyahBott.AtiyahBottFormula","page":"AtiyahBott.jl","title":"AtiyahBott.AtiyahBottFormula","text":"AtiyahBottFormula(n, d, m, P, do_check, show_bar, down_col)\n\nApply the Atiyah-Bott residue formula to the class P, in the moduli space of rational marked stable maps to the projective space of dimension n of degree d with m marks.\n\nArguments\n\nn::Int64: the dimension of the projective space, it must be between 1 and 254.\nd::Int64: the degree of the stable maps, it must be between 1 and 13.\nm::Int64: the number of marks.\nP: the equivariant class.\ndo_check::Bool: if true, checks if P is a well defined zero cycle, and stops the computation if this is not true. If false, the computation may have an unexpected behaviour. By default is true.\nshow_bar::Bool: hide the progress bar if and only if this condition is false. By default is true.\ndown_col::Bool: check if all colorations needed in the computation are in the folder Data, and download them otherwise. Once downloaded, those files can be used for future computations. By default is true.\n\nThe general construction of P is the following:\n\njulia> P = (g,c,w,s,m) ->\n\nAfter ->, one has to write an expression in the equivariant classes. All such equivariant classes are functions starting with (g,c,w,s) or (g,c,w,s,m). At the end, they can have more arguments. The expression is a polynomial combination of the equivariant classes. We compute the degree of P by\n\njulia> AtiyahBottFormula(n,d,m,P);\n\nExample\n\njulia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);\njulia> AtiyahBottFormula(3,1,0,P);\nWarning: the class is not a 0-cycle.\njulia> AtiyahBottFormula(4,1,0,P);\nResult: 2875//1\njulia> AtiyahBottFormula(4,1,0,P,false);             #same as before, but without the preliminary check on `P`\njulia> AtiyahBottFormula(4,1,0,P,false,false);       #same as before, but without showing the progress bar\njulia> AtiyahBottFormula(4,1,0,P,false,false,false); #same as before, but without checking for the colorations\n\nThe function returns an array of the same dimension of P (non-vectorized classes are assumed as 1-dimensional arrays). The Julia notation for accessing to array is name_of_array[i] where i is an index starting from 1.\n\nExample\n\njulia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)*Hypersurface(g,c,w,s,3);\njulia> x = AtiyahBottFormula(3,2,0,P)[1];\nResult: 81//1\njulia> x\n81//1\n\nThe class P supports parameters.\n\njulia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,3)*(Incidency(g,c,w,s,2)//3)^(d-1);\njulia> d = 2;\njulia> AtiyahBottFormula(3,d,0,P);\nResult: 27//1\njulia> d = 3;\njulia> AtiyahBottFormula(3,d,0,P);\nResult: 84//1\n\nMore examples are available in the support of the equivariant classes. It is enough to type ? and then the name of the class. Currently, the supported classes are:\n\nO1_i         (Euler class of mathrmev_i^*mathcalO_mathbbP^n(1))\nO1           (product of all O1_i)\nPsi          (cycle of psi-classes)\nJet          (Euler class of the jet bundle J^p of mathrmev^*mathcalO_mathbbP^n(q))\nHypersurface (Euler class of the direct image of mathrmev^*mathcalO_mathbbP^n(b))\nIncidency    (cycle parameterizing curves meeting a linear subspace)\nContact      (cycle parameterizing contact curves)\nR1           (first derived functor of direct image of mathrmev^*mathcalO_mathbbP^n(-k))\n\nTo add more classes, please contact the authors.\n\n\n\n\n\n","category":"function"},{"location":"#Equivariant-Classes","page":"AtiyahBott.jl","title":"Equivariant Classes","text":"","category":"section"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"Here we list all equivariant classes currently supported by the package.","category":"page"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"Hypersurface\r\nO1_i\r\nO1\r\nIncidency\r\nPsi\r\nJet\r\nContact\r\nR1\r\nEuler_inv","category":"page"},{"location":"#AtiyahBott.Hypersurface","page":"AtiyahBott.jl","title":"AtiyahBott.Hypersurface","text":"Hypersurface(g, c, w, s, b)\n\nEquivariant class of the Euler class of the bundle equal to the direct image under the forgetful map of mathrmev^*mathcalO_mathbbP^n(b). It parameterizes curves contained in a hypersurface of degree b.\n\nArguments\n\ng::SimpleGraph: the graph.\nc::Vector{UInt8}: the coloration.\nw::Vector{Int64}: the weights.\ns::Rational{BigInt}: the scalars.\nb::Int64: the degrees of the hypersurface. Alternatively, it can be an array of integers, meaning the multiplication of the equivariant class defined by each element of the array.\n\nExample\n\nThe following Gromov-Witten invariants of Calabi-Yau threefolds\n\nbeginaligned\nint_overlineM_00(mathbbP^41)mathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^4(5))) = 2875 \nint_overlineM_00(mathbbP^52)mathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^5(3))^oplus 2) = frac4235498 \nint_overlineM_00(mathbbP^53)mathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^5(4)))cdotmathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^5(2))) = frac42269081627 \nint_overlineM_00(mathbbP^74)mathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^7(2)))^4 = 25705160 \nendaligned\n\ncan be computed as\n\njulia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);\njulia> AtiyahBottFormula(4,1,0,P);\nResult: 2875//1\njulia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,[3,3]);\njulia> AtiyahBottFormula(5,2,0,P);\nResult: 423549//8\njulia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,4)*Hypersurface(g,c,w,s,2);\njulia> AtiyahBottFormula(5,3,0,P);\nResult: 422690816//27\njulia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,2)^4;\njulia> AtiyahBottFormula(7,4,0,P);\nResult: 25705160//1\n\nwarning: Attention!\nThe program will stop if b is not positive.\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.O1_i","page":"AtiyahBott.jl","title":"AtiyahBott.O1_i","text":"O1_i(g, c, w, s, m, i)\n\nEquivariant class of the pull-back of mathcalO_mathbbP^n(1) with respect to the i-th evaluation map.\n\nArguments\n\ng::SimpleGraph: the graph.\nc::Vector{UInt8}: the coloration.\nw::Vector{Int64}: the weights.\ns::Rational{BigInt}: the scalars.\nm::Marks: the marks.\ni::Int64: the evaluation map.\n\nExample\n\nThe following Gromov-Witten invariants\n\nbeginaligned\nint_overlineM_02(mathbbP^21)mathrmev_1^*mathcalO_mathbbP^2(1)cdotmathrmev_2^*mathcalO_mathbbP^2(1) = 1 \nint_overlineM_01(mathbbP^31)mathrmev_1^*mathcalO_mathbbP^2(1)^2cdotmathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^3(2))) = 4\nendaligned\n\ncan be computed as\n\njulia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)*O1_i(g,c,w,s,m,2);\njulia> AtiyahBottFormula(2,1,2,P);\nResult: 1//1\njulia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)^2*Hypersurface(g,c,w,s,2);\njulia> AtiyahBottFormula(3,1,1,P);\nResult: 4//1\n\nwarning: Attention!\nThe program will stop if i is not between 1 and the number of marks.\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.O1","page":"AtiyahBott.jl","title":"AtiyahBott.O1","text":"O1(g, c, w, s, m)\n\nEquivariant class of the pull-back of mathcalO_mathbbP^n(1) with respect to the product of all evaluation maps.\n\nArguments\n\ng::SimpleGraph: the graph.\nc::Vector{UInt8}: the coloration.\nw::Vector{Int64}: the weights.\ns::Rational{BigInt}: the scalars.\nm::Marks: the marks.\n\nThis function is equivalent to the product of the function O1_i(g,c,w,s,m,i) where i runs from 1 to the number of marks.\n\nExample\n\nThe following Gromov-Witten invariants\n\nbeginaligned\nint_overlineM_08(mathbbP^23)prod_i=1^8mathrmev_i^*mathcalO_mathbbP^2(1) = 12 \nint_overlineM_01(mathbbP^32)mathrmev_1^*mathcalO_mathbbP^3(1)^2cdotmathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^3(3))) = 81 \nendaligned\n\ncan be computed as\n\njulia> P = (g,c,w,s,m) -> O1(g,c,w,s,m);\njulia> AtiyahBottFormula(2,3,8,P);\nResult: 12//1\njulia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Hypersurface(g,c,w,s,3);\njulia> AtiyahBottFormula(3,2,1,P);\nResult: 81//1\n\nIn order to remove O1_i(g,c,w,s,m,i) for some i, it is enough to divide by that function.\n\nExample\n\njulia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)//O1_i(g,c,w,s,m,1);\n\nHere P is the product of all O1_i(g,c,w,s,m,i) where i runs from 2 to m.\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.Incidency","page":"AtiyahBott.jl","title":"AtiyahBott.Incidency","text":"Incidency(g, c, w, s, r)\n\nEquivariant class of the cycle parameterizing curves meeting a linear subspace of codimension r.\n\nArguments\n\ng::SimpleGraph: the graph.\nc::Vector{UInt8}: the coloration.\nw::Vector{Int64}: the weights.\ns::Rational{BigInt}: the scalars.\nr::Int64: the codimension of the subvariety. Alternatively, it can be an array of integers, meaning the multiplication of the equivariant class defined by each element of the array.\n\nExample\n\nThe following Gromov-Witten invariants\n\nbeginaligned\nint_overlineM_00(mathbbP^31)delta_*(mathrmev^*mathcalO_mathbbP^3(1)^3)^2 = 1 \nint_overlineM_00(mathbbP^31)delta_*(mathrmev^*mathcalO_mathbbP^3(1)^2)^2cdot delta_*(mathrmev^*mathcalO_mathbbP^3(1)^3) = 1 \nint_overlineM_00(mathbbP^33)delta_*(mathrmev^*mathcalO_mathbbP^3(1)^2)^2cdot mathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^3(3))) = 756 \nendaligned\n\ncan be computed as\n\njulia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,3)^2;\njulia> AtiyahBottFormula(3,1,0,P);\nResult: 1//1\njulia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,[2,2,3]);\njulia> AtiyahBottFormula(3,1,0,P);\nResult: 1//1\njulia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,[2,2])*Hypersurface(g,c,w,s,3);\njulia> AtiyahBottFormula(3,3,0,P);\nResult: 756//1\n\nwarning: Attention!\nThe program will stop if r is not positive.\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.Psi","page":"AtiyahBott.jl","title":"AtiyahBott.Psi","text":"Psi(g, c, w, s, m, a)\n\nEquivariant class of the cycle of psi-classes.\n\nArguments\n\ng::SimpleGraph: the graph.\nc::Vector{UInt8}: the coloration.\nw::Vector{Int64}: the weights.\ns::Rational{BigInt}: the scalars.\nm::Marks: the marks.\na::Vector{Int64}: the vector of the exponents of the psi classes. It is ordered, meaning that the first element is the exponent of psi_1, the second is the exponent of psi_2, and so on.\n\nnote: Note\nThe size of a must be at most m. If it is smaller, missing exponents will be considered as zeros. If a is a number, it will be considered as the exponent of psi_1.\n\nwarning: Attention!\nThe program will stop if we have one of the following conditions:the size of a is bigger than m,\na contains a negative number.\n\nExample\n\nThe following Gromov-Witten invariants\n\nbeginaligned\nint_overlineM_02(mathbbP^62)mathrmev_1^*mathcalO_mathbbP^6(1)^5cdotmathrmev_2^*mathcalO_mathbbP^6(1)^2cdotmathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^6(5)))cdotpsi_1psi_2^0 = 495000 \nint_overlineM_02(mathbbP^102)mathrmev_1^*mathcalO_mathbbP^10(1)^8cdotmathrmev_2^*mathcalO_mathbbP^10(1)^6cdotmathrmc_top(delta_*(mathrmev^*mathcalO_mathbbP^10(7)))cdotpsi_1^2 = 71804533752 \nint_overlineM_01(mathbbP^22)mathrmev_1^*mathcalO_mathbbP^2(1)^2cdotpsi_1^4 = frac18 \nint_overlineM_01(mathbbP^22)delta_*(mathrmev^*mathcalO_mathbbP^2(1)^2)^4cdotmathrmev_1^*mathcalO_mathbbP^2(1)cdot(mathrmev_1^*mathcalO_mathbbP^2(1)+psi_1) = 2 \nint_overlineM_01(mathbbP^22)mathrmev_1^*mathcalO_mathbbP^2(1)^2cdot(psi_1^3cdotmathrmev_1^*mathcalO_mathbbP^2(1)+psi_1^2cdotmathrmev_1^*mathcalO_mathbbP^2(1)^2) = frac18 \nint_overlineM_01(mathbbP^32)mathrmev_1^*mathcalO_mathbbP^2(1)cdot(psi_1^7cdotmathrmev_1^*mathcalO_mathbbP^2(1)+psi_1^6cdotmathrmev_1^*mathcalO_mathbbP^2(1)^2) = -frac516 \nendaligned\n\ncan be computed as\n\njulia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)^5*O1_i(g,c,w,s,m,2)^2*Hypersurface(g,c,w,s,5)*Psi(g,c,w,s,m,[1,0]);\njulia> AtiyahBottFormula(6,2,2,P);\nResult: 495000//1\njulia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)^8*O1_i(g,c,w,s,m,2)^6*Hypersurface(g,c,w,s,7)*Psi(g,c,w,s,m,2);\njulia> AtiyahBottFormula(10,2,2,P);\nResult: 71804533752//1\njulia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Psi(g,c,w,s,m,4);\njulia> AtiyahBottFormula(2,2,1,P);\nResult: 1//8\njulia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^4*O1_i(g,c,w,s,m,1)*(O1_i(g,c,w,s,m,1) + Psi(g,c,w,s,m,1))\njulia> AtiyahBottFormula(2,2,1,P); #number of plane conics through four points and tangent to a line\nResult: 2\njulia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*(Psi(g,c,w,s,m,3)*O1(g,c,w,s,m)+Psi(g,c,w,s,m,2)*O1(g,c,w,s,m)^2);\njulia> AtiyahBottFormula(2,2,1,P);\nResult: 1//8\njulia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)*(Psi(g,c,w,s,m,7)*O1(g,c,w,s,m)+Psi(g,c,w,s,m,6)*O1(g,c,w,s,m)^2);\njulia> AtiyahBottFormula(3,2,1,P);\nResult: -5//16\n\nwarning: Psi is singleton!\nPsi cannot be multiplied by itself.julia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Psi(g,c,w,s,m,1)^4;                  #this is **wrong**\njulia> AtiyahBottFormula(2,2,1,P);\nWarning: more instances of Psi has been found. Type:\njulia> ?Psi\nfor support.\njulia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Psi(g,c,w,s,m,3)*Psi(g,c,w,s,m,1);   #this is **wrong**\njulia> AtiyahBottFormula(2,2,1,P);\nWarning: more instances of Psi has been found. Type:\njulia> ?Psi\nfor support.\njulia> P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2*Psi(g,c,w,s,m,4);\njulia> AtiyahBottFormula(2,2,1,P);\nResult: 1//8\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.Jet","page":"AtiyahBott.jl","title":"AtiyahBott.Jet","text":"Jet(g, c, w, s, m, p, q)\n\nEquivariant class of the jet bundle J^p of the pull back of mathcalO_mathbbP^n(q) with respect to the first psi-class.\n\nArguments\n\ng::SimpleGraph: the graph.\nc::Vector{UInt8}: the coloration.\nw::Vector{Int64}: the weights.\ns::Rational{BigInt}: the scalars.\nm::Marks: the marks.\np::Int64: the exponent of the Jet bundle. In particular, it is a bundle of rank p+1.\nq::Int64: the degree of the line bundle that is pulled back.\n\nnote: Note\nIn order to define this bundle, the number of marks must be at least 1. You cannot multiply this bundle by the class Psi(g,c,w,s,m,a).\n\nExample\n\nbeginaligned\nint_overlineM_01(mathbbP^22)delta_*(mathrmev^*mathcalO_mathbbP^2(1)^2)^4cdotmathrmc_top(J^1(mathrmev_1^*mathcalO_mathbbP^2(1))) = 2 \nint_overlineM_01(mathbbP^22)delta_*(mathrmev^*mathcalO_mathbbP^2(1)^2)^4cdot(mathrmc_top(J^1(mathrmev_1^*mathcalO_mathbbP^2(1)))+mathrmev_1^*mathcalO_mathbbP^2(1)^2) = 3 \nint_overlineM_01(mathbbP^3d)fracmathrmev_1^*mathcalO_mathbbP^3(1)^2kcdotmathrmc_top(J^4d-2(mathrmev_1^*mathcalO_mathbbP^3(k))) = frac(4d-2)(d)^4 \nendaligned\n\ncan be computed as\n\njulia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^4*Jet(g,c,w,s,m,1,1);\njulia> AtiyahBottFormula(2,2,1,P);\nResult: 2//1\njulia> P = (g,c,w,s,m) -> Incidency(g,c,w,s,2)^4*(Jet(g,c,w,s,m,1,1)+O1(g,c,w,s,m)^2);\njulia> AtiyahBottFormula(2,2,1,P);\nResult: 3//1\njulia> P = (g,c,w,s,m) -> (O1(g,c,w,s,m)^2)//k*Jet(g,c,w,s,m,4*d-2,k);\njulia> d=1;k=1;AtiyahBottFormula(3,d,1,P);   #The value of this integral does not depend on k, only on d\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.Contact","page":"AtiyahBott.jl","title":"AtiyahBott.Contact","text":"Contact(g, c, w, s)\n\nEquivariant class of the Euler class of the bundle equal to the direct image under the forgetful map of mathrmev^*mathcalO_mathbbP^n(2) tensor the dualizing sheaf of the forgetful map. It parameterizes contact curves in an odd dimensional projective space.\n\nArguments\n\ng::SimpleGraph: the graph.\nc::Vector{UInt8}: the coloration.\nw::Vector{Int64}: the weights.\ns::Rational{BigInt}: the scalars.\n\nExample\n\nbeginaligned\nint_overlineM_02(mathbbP^31)mathrmev_1^*mathcalO_mathbbP^3(1)^2cdotmathrmev_2^*mathcalO_mathbbP^3(1)^3cdotmathrmc_top(delta_*(omega_deltaotimesmathrmev^*mathcalO_mathbbP^3(2))) = 1 \nendaligned\n\ncan be computed as\n\njulia> P = (g,c,w,s,m) -> O1_i(g,c,w,s,m,1)^2*O1_i(g,c,w,s,m,2)^3*Contact(g,c,w,s);\njulia> AtiyahBottFormula(3,1,2,P);\nResult: 1//1\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.R1","page":"AtiyahBott.jl","title":"AtiyahBott.R1","text":"R1(g, c, w, s, k)\n\nThe equivariant class of the first derived functor of the pull-back of mathcalO_mathbbP^n(-k).\n\nArguments\n\ng::SimpleGraph: the graph.\nc::Vector{UInt8}: the coloration.\nw::Vector{Int64}: the weights.\ns::Rational{BigInt}: the scalars.\nk::Int64: a positive integer.\n\nExample\n\nbeginaligned\nint_overlineM_00(mathbbP^1d)mathrmc_top(R^1delta_*(mathrmev^*mathcalO_mathbbP^3(-1)))^2 = frac1d^3 \nendaligned\n\ncan be computed as\n\njulia> d = 1; #for other values of d, change this line\njulia> P = (g,c,w,s,m) -> R1(g,c,w,s,1)^2;\njulia> AtiyahBottFormula(1,d,0,P);\nResult: 1//1\n\nwarning: Attention!\nThe program will stop if k is not positive.\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.Euler_inv","page":"AtiyahBott.jl","title":"AtiyahBott.Euler_inv","text":"Euler_inv(g, c, w, s, m)\n\nThe inverse of the (equivariant) Euler class of the normal bundle. This function is invoked automatically.\n\nArguments\n\ng::SimpleGraph: the graph.\nc::Vector{UInt8}: the coloration.\nw::Vector{Int64}: the weights.\ns::Rational{BigInt}: the scalars.\nm::Marks: the marks.\n\n\n\n\n\n","category":"function"},{"location":"#Other-Functions","page":"AtiyahBott.jl","title":"Other Functions","text":"","category":"section"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"dim_M\r\ncodim\r\nis_zero_cycle\r\ncheck_Data\r\nfill_Data\r\nfree_Data\r\nAtiyahBottFormulaForGraph","category":"page"},{"location":"#AtiyahBott.dim_M","page":"AtiyahBott.jl","title":"AtiyahBott.dim_M","text":"dim_M(n, d, m)\n\nThe dimension of the moduli space of stable rational map to the projective space of dimension n, of degree d with m marks.\n\nArguments\n\nn::Int64: the dimension of the projective space.\nd::Int64: the degree of the stable maps.\nm::Int64: the number of marks.\n\nExample\n\njulia> dim_M(2,2,5)\n10\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.codim","page":"AtiyahBott.jl","title":"AtiyahBott.codim","text":"codim(n, d, m, P)\n\nThe codimension of the equivariant class P.\n\nArguments\n\nn::Int64: the dimension of the projective space.\nd::Int64: the degree of the stable maps.\nm::Int64: the number of marks.\nP: the equivariant class.\n\nExample\n\njulia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);\njulia> codim(4,1,0,P)\n6\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.is_zero_cycle","page":"AtiyahBott.jl","title":"AtiyahBott.is_zero_cycle","text":"is_zero_cycle(n, d, m, P)\n\nReturn true if the equivariant class P is a 0-cycle in the moduli space, false otherwise.\n\nArguments\n\nn::Int64: the dimension of the projective space.\ndeg::Int64: the degree of the stable maps.\nn_marks::Int64: the number of marks.\nP: the equivariant class.\n\nExample\n\njulia> P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5);\njulia> is_zero_cycle(4,1,0,P)\ntrue\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.check_Data","page":"AtiyahBott.jl","title":"AtiyahBott.check_Data","text":"check_Data()\n\nList of all files containing the colorations in the folder Data.\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.fill_Data","page":"AtiyahBott.jl","title":"AtiyahBott.fill_Data","text":"fill_Data(n, d)\n\nDownload from internet all colorations used for computations in the moduli space with dimension n and degree d. Return true if there is no need to download any file or if all downloads go well. Otherwise, return false.\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.free_Data","page":"AtiyahBott.jl","title":"AtiyahBott.free_Data","text":"free_Data()\n\nDelete the folder Data.\n\n\n\n\n\n","category":"function"},{"location":"#AtiyahBott.AtiyahBottFormulaForGraph","page":"AtiyahBott.jl","title":"AtiyahBott.AtiyahBottFormulaForGraph","text":"AtiyahBottFormulaForGraph(g, pruf_str, aut, n, deg, n_marks, P, s)\n\nApply the Atiyah-Bott formula to all colorations of a specific graph. It is useful for splitting the computation in multiple parts, to be computed in single threads.\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"AtiyahBott.jl","title":"Index","text":"","category":"section"},{"location":"","page":"AtiyahBott.jl","title":"AtiyahBott.jl","text":"","category":"page"}]
}
