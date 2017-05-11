module PARI

export @gp_str, pari, getavma, PariInt, Gen, GEN

include("./Utils.jl")
include("./Lib.jl")
include("./GP.jl")

using .Utils
using .GP
using .Lib
using .Lib: CulongMax, ClongMax, at0, at0!, at1, at1!





function pari_print(io::IO, a::Gen)
    cstr = GENtostr(a.ptr)
    print(io, unsafe_string(cstr))
    pari_free(cstr)
end

Base.show(io::IO, ::MIME{Symbol("text/plain")}, g::Gen) = pari_print(io, g)
Base.show(io::IO, g::Gen) = pari_print(io, g)
Base.showcompact(io::IO, g::Gen) = show(io, g)



"Make an unitialized INT encodable with `s` Limbs,
initializing the 2-words header"
function make_int(s::Integer, signe::Union{Int64,Int32})
    s >= 0 || throw(DomainError())
    ss = s+2
    x = Gen(ss)
    at0!(x.ptr, evaltyp(t_INT) | evallg(ss))
    at1!(x.ptr, evalsigne(signe) | evallgefint(ss))
    return x
end

function pari(a::BigInt)
    s = abs(a.size)
    ss = s+2
    x = make_int(s, sign(a.size))
    for i in 1:s
#       unsafe_store!(int_W_lg(x.ptr, i-1, ss),
#                     unsafe_load(a.d, i) % long)
        unsafe_store!(x.ptr,
                      unsafe_load(a.d, i) % long,
                      #ss-i+1) # WARNING: valid only for native kernel!
                      2+i)
    end
    return x
end

@assert UInt === ulong && Int === long === Clong

function pari(a::CulongMax)
    ss = sign(a) % Clong
    x = make_int(ss, ss)
    (ss % Bool) && unsafe_store!(x.ptr, a % Clong, 3)
    x
end

function pari(a::ClongMax)
    x = pari(abs(a) % Culong)
    a < 0 && setsigne(x.ptr, sign(a))
    x
end

pari(a::Integer) = pari(big(a))


import Base: +, *, -, div, rem

function +(x::Gen, y::Gen)
    typ(x) > t_REAL || typ(y) > t_REAL && throw(DomainError())
    mpadd(x, y)
end

struct PariInt <: Integer
    i::Gen

    function PariInt(x::Gen)
        typ(x) == t_INT || throw(DomainError())
        new(x)
    end
end

PariInt(x::Integer) = PariInt(pari(x))
PariInt(x::PariInt) = x

# TODO: assuming GMP
function Base.convert(::Type{long}, x::Gen)
    typ(x) == t_INT || throw(InexactError())
    lg(x) < 3 && return 0
    lgefint(x.ptr) > 3 && throw(InexactError())
    y = unsafe_load(x.ptr, 3)
    signe(x.ptr) > 0 ? y : -y
end

Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::PariInt) = pari_print(io, x.i)



+(x::PariInt, y::PariInt) = PariInt(x.i + y.i)

factor(x::Integer) = let f = convert(Array, GP.factor(PariInt(x).i)),
    d = Dict{Gen, Int}()
    n = size(f, 1)
    for i = 1:n
        d[f[i, 1]] = f[i, 2]
    end
    d
end

# WARNING: this does NOT take ownership of elements
# We should do a copy to avoid GC making those references unreachable
function Base.convert(::Type{Array}, x::Gen)
    n = lg(x)-1
    if typ(x) == t_COL
        a = Vector{Gen}(n)
        for i = i:n
            a[i] = x[i]
        end
    elseif typ(x) == t_VEC
        a = Array{Gen}(1, n)
        for i = i:n
            a[i] = x[i]
        end
    elseif typ(x) == t_MAT
        n < 1 && return Array{Gen}(0, 0)
        m = lg(x[1])-1
        a = Array{Gen}(m, n)
        for j = 1:n,
            i = 1:m
            a[i, j] = x[i, j]
        end
    end
    return a
end

end # module pari
