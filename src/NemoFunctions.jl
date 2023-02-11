###############################################################################
# 
# low-level in-place arithmetic operators
# should probably be added into Nemo upstream
# 
#=
function set!(x::fmpz, y::Int)
    ccall((:fmpz_set_si, Nemo.libflint), Nothing, (Ref{fmpz}, Int), x, y)
    return x
end

function set!(x::fmpq, y::fmpq)
    ccall((:fmpz_set, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}), x, y)
    return x
end

function set!(x::fmpq, n::Int, d::Int)
        ccall((:fmpq_set_si, Nemo.libflint), Nothing, (Ref{fmpq}, Int, UInt), x, n, UInt(d))
    return x
end

function set!(x::fmpq, n::fmpz, d::fmpz)
    ccall((:fmpq_set_fmpz_frac, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpz}, Ref{fmpz}), x, n, d)
    return x
end

function Nemo.mul!(z::fmpq, x::fmpq, y::fmpz)
    ccall((:fmpq_mul_fmpz, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}, Ref{fmpz}), z, x, y)
    return z
end

function Nemo.mul!(z::fmpq, x::fmpq, y::Int)
    ccall((:fmpq_mul_si, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}, Int), z, x, y)
    return z
end

function Nemo.sub!(z::fmpq, x::fmpq, y::Int)
    ccall((:fmpq_sub_si, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}, Int), z, x, y)
    return z
end

function Nemo.inv!(x::fmpq, y::fmpq)
    ccall((:fmpq_inv, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}), x, y)
    return x
end

function pow!(c::fmpq, a::fmpq, b::Int)
    iszero(a) && b < 0 && throw(DivideError())
    ccall((:fmpq_pow_si, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}, Int), c, a, b)
    return c
end

function div!(z::fmpq, x::fmpq, y::fmpq)
    iszero(y) && throw(DivideError())
    ccall((:fmpq_div, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}), z, x, y)
    return z
end
  
muleq!(x::fmpq, y::fmpq) = mul!(x, x, y)
muleq!(x::fmpq, y::fmpz) = mul!(x, x, y)
muleq!(x::fmpq, y::Int) = mul!(x, x, y)
subeq!(x::fmpq, y::Int) = Nemo.sub!(x, x, y)
poweq!(x::fmpq, y::Int) = pow!(x, x, y)
diveq!(x::fmpq, y::fmpq) = div!(x, x, y)
Nemo.inv!(x::fmpq) = Nemo.inv!(x, x)
=#
"""
    pow_eq!(a, b)
Equivalent to a ^= b
"""
function pow_eq!(a::fmpq, b::Int)::Nothing
    #iszero(a) && b < 0 && throw(DivideError())
    ccall((:fmpq_pow_si, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}, Int), a, a, b)
end

"""
    mul_eq!(a, b)
Equivalent to a *= b
"""
# function mul_eq!(a::fmpq, b::fmpz)::Nothing
#     ccall((:fmpq_mul, Nemo.libflint), Nothing,(Ref{fmpq}, Ref{fmpq}, Ref{fmpz}), a, a, b)
# end
function mul_eq!(a::fmpq, b::fmpq)::Nothing
    ccall((:fmpq_mul, Nemo.libflint), Nothing,(Ref{fmpq}, Ref{fmpq}, Ref{fmpq}), a, a, b)
end
function mul_eq!(a::fmpq, b::Int)::Nothing
    ccall((:fmpq_mul_si, Nemo.libflint), Nothing,(Ref{fmpq}, Ref{fmpq}, Int), a, a, b)
end

"""
    add_eq!(a, b)
Equivalent to a += b
"""
function add_eq!(a::fmpq, b::fmpq)::Nothing
    ccall((:fmpq_add, Nemo.libflint), Nothing,(Ref{fmpq}, Ref{fmpq}, Ref{fmpq}), a, a, b)
end
function add_eq!(a::Vector{fmpq}, b::Vector{fmpq})::Nothing
    for i in 1:length(a)
        ccall((:fmpq_add, Nemo.libflint), Nothing,(Ref{fmpq}, Ref{fmpq}, Ref{fmpq}), a[i], a[i], b[i])
    end
end

"""
    div_eq!(a, b)
Equivalent to a /= b
"""
function div_eq!(a::fmpq, b::fmpq)::Nothing
    ccall((:fmpq_div, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}), a, a, b)
end
function div_eq!(a::fmpq, b::Int)::Nothing
    ccall((:fmpq_div_fmpz, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}, Ref{fmpz}), a, a, fmpz(b))
end

"""
    eq!(a, b)
Equivalent to a = b
"""
function eq!(a::fmpq, b::fmpq)::Nothing
    ccall((:fmpq_set, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}), a, b)
end

"""
    sub!(a, b)
Equivalent to a -= b
"""
function sub!(a::fmpq, b::fmpq)::Nothing
    ccall((:fmpq_sub, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}, Ref{fmpq}), a, a, b)
end

"""
    neg!(a, b)
Equivalent to a *= -1
"""
function neg!(a::fmpq)::Nothing
    ccall((:fmpq_neg, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}), a, a)
end

"""
    inv!(a, b)
Equivalent to a /= 1
"""
function inv!(a::fmpq)::Nothing
    ccall((:fmpq_inv, Nemo.libflint), Nothing, (Ref{fmpq}, Ref{fmpq}), a, a)
end