# AtiyahBott.jl


```@contents
```

```@docs
AtiyahBott
```
In order to install the package, type:
```julia-repl
julia> using Pkg
julia> Pkg.add("AtiyahBott")
```
To load the package, type:
```julia-repl
julia> using AtiyahBott
```
## The function AtiyahBottFormula
This is the main function of the package.
```@docs
AtiyahBottFormula(n::Int64, deg::Int64, n_marks::Int64, P, do_check::Bool = true, show_bar::Bool = true, down_col::Bool = true)
```

## Equivariant Classes
Here we list all equivariant classes currently supported by the package.
```@docs
Hypersurface
O1_i
O1
Incidency
Psi
Jet
Contact
R1
Euler_inv
```

## Other Functions
```@docs
dim_M
codim
is_zero_cycle
check_Data
fill_Data
free_Data
```

## Index

```@index
```
