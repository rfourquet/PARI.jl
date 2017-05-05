module PARI

export @pari_str, pari, getavma, PariInt, Gen, GEN

using Base: unsafe_convert
using Base.GMP: CulongMax, ClongMax

# macros
const LONG_IS_64_BIT = sizeof(Clong)*8 == 64

@assert Culonglong === Culong
const pari_ulong = Culong
const ulong = Culong
const long = Clong

const IUL = one(ulong) # == 1UL in C

if LONG_IS_64_BIT
    const BITS_IN_LONG = 64
    const BYTES_IN_LONG = 8 # not in PARI
    const TWOPOTBITS_IN_LONG = 6
    const LONG_MAX = Clong(9223372036854775807)
    SMALL_ULONG(p) = (p % ulong) <= (3037000493 % ulong)
else
    const BITS_IN_LONG = 32
    const BYTES_IN_LONG = 4 # not in PARI
    const TWOPOTBITS_IN_LONG = 5
    const LONG_MAX = Clong(2147483647)
    SMALL_ULONG(p) = (p % ulong) <= (46337 % ulong)
end
const ULONG_MAX = ~(zero(ulong))

const DEFAULTPREC    = (2 + 8  ÷ sizeof(long))
const MEDDEFAULTPREC = (2 + 16 ÷ sizeof(long))
const BIGDEFAULTPREC = (2 + 24 ÷ sizeof(long))
const LOWDEFAULTPREC =  3
const EXTRAPRECWORD  =  1
const HIGHBIT = IUL << (BITS_IN_LONG-1)
const BITS_IN_HALFULONG = BITS_IN_LONG >> 1

const LOWMASK = IUL << BITS_IN_HALFULONG - IUL
const HIGHMASK = ~LOWMASK

HIGHWORD(a) = a >> BITS_IN_HALFULONG
LOWWORD(a) = a & LOWMASK

const TYPnumBITS  = 7
const SIGNnumBITS = 2

const VARNnumBITS = LONG_IS_64_BIT ? 16 : 14

const   LGnumBITS = BITS_IN_LONG - 1 - TYPnumBITS
const VALPnumBITS = BITS_IN_LONG - SIGNnumBITS - VARNnumBITS
const EXPOnumBITS = BITS_IN_LONG - SIGNnumBITS
const PRECPSHIFT = VALPnumBITS
const  VARNSHIFT = VALPnumBITS
const   TYPSHIFT = BITS_IN_LONG - TYPnumBITS
const  SIGNSHIFT = BITS_IN_LONG - SIGNnumBITS

const EXPOBITS    = (IUL << EXPOnumBITS) - 1
const SIGNBITS    = ~((IUL << SIGNSHIFT) - 1)
const  TYPBITS    = ~((IUL <<  TYPSHIFT) - 1)
const VALPBITS    = (IUL << VALPnumBITS) - 1
const PRECPBITS   = ~VALPBITS
const LGBITS      = (IUL << LGnumBITS) - 1
const MAXVARN     = (IUL << VARNnumBITS) - 1
const VARNBITS    = MAXVARN << VARNSHIFT
const NO_VARIABLE = 2147483647 % long

const HIGHEXPOBIT = IUL << (EXPOnumBITS-1)
const HIGHVALPBIT = IUL << (VALPnumBITS-1)
const CLONEBIT    = IUL << LGnumBITS

evaltyp(x)      = (x % ulong) << TYPSHIFT
evalvarn(x)     = (x % ulong) << VARNSHIFT
evalsigne(x)    = ((x % long) << SIGNSHIFT) % ulong
_evalexpo(x)    = HIGHEXPOBIT + x
_evalvalp(x)    = HIGHVALPBIT + x
_evalprecp(x)   = (x % long) << PRECPSHIFT
evallgefint(x)  = x
evallgeflist(x) = x
_evallg(x)      = x

typ(x) = ((unsafe_load(x) % ulong) >> TYPSHIFT) % long

at0(x) = unsafe_load(x)
at(x) = unsafe_load(x)
at1(x) = unsafe_load(x, 2)
at0!(x, a) = unsafe_store!(x, a)
at!(x, a) = unsafe_store!(x, a)
at1!(x, a) = unsafe_store!(x, a, 2)
ul_(x) = unsafe_convert(Ptr{ulong}, x)

settyp(x, s) = let xx = ul_(x)
    at0!(xx, (at0(xx) & (~TYPBITS)) | evaltyp(s))
end

isclone(x) = at0(ul_(x)) & CLONEBIT
setisclone(x) = at!(ul_(x), at(ul_(x)) | CLONEBIT)
unsetisclone(x) = at!(ul_(x), at(ul_(x)) & ~CLONEBIT)

lg(x) = ((at(x) % ulong) & LGBITS) % long
setlg(x, s) = at!(ul_(x), (at(ul_(x)) & (~LGBITS)) | evallg(s))

signe(x) = (at1(x) % long) >> SIGNSHIFT
setsigne(x, s) = at1!(ul_(x),
                      (at1(ul_(x)) & (~SIGNBITS)) | (evalsigne(s) % ulong))

lgefint(x) = ((at1(x) % ulong) & LGBITS) % long
setlgefint(x, s) = at1!(ul_(x),
                      (at1(ul_(x)) & (~LGBITS)) | (evallgefint(s) % ulong))


# Warning: for native kernel
int_W(x, l) = x+(lgefint(x)-1-l)*BYTES_IN_LONG
int_W_lg(x, l, lx) = x+(lx-1-l)*BYTES_IN_LONG



const t_INT    =  1
const t_REAL   =  2
const t_INTMOD =  3
const t_FRAC   =  4
const t_FFELT  =  5
const t_COMPLEX=  6
const t_PADIC  =  7
const t_QUAD   =  8
const t_POLMOD =  9
const t_POL    =  10
const t_SER    =  11
const t_RFRAC  =  13
const t_QFR    =  15
const t_QFI    =  16
const t_VEC    =  17
const t_COL    =  18
const t_MAT    =  19
const t_LIST   =  20
const t_STR    =  21
const t_VECSMALL= 22
const t_CLOSURE = 23
const t_ERROR   = 24


@enum(Types,
      INT=1,
      REAL,
      INTMOD,
      FRAC,
      FFELT,
      COMPLEX,
      PADIC,
      QUAD,
      POLMOD,
      POL,
      SER,
      RFRAC=13,
      QFR=15,
      QFI,
      VEC,
      COL,
      MAT,
      LIST,
      STR,
      VECSMALL,
      CLOSURE,
      ERROR)

# end macros

const libpari = :libpari

push!(Libdl.DL_LOAD_PATH, "/usr/local/lib")
Libdl.dlopen(libpari)

pari_init(size::Integer, maxprime::Integer) =
    ccall((:pari_init, libpari), Void, (Csize_t, Culong), size, maxprime)

pari_close() = ccall((:pari_close, libpari), Void, ())

pari_malloc(size::Integer) =
    ccall((:pari_malloc, libpari), Ptr{Void}, (Csize_t,), size)

pari_free(ptr::Ptr) =
    ccall((:pari_free, libpari), Void, (Ptr{Void},), ptr)


evallg(x) = ccall((:evallg, libpari), long, (long,), x)




const GEN = Ptr{Clong}

mutable struct Gen # mutable for allowing finalizers
    ptr::GEN

    function Gen(ptr::GEN, clone::Bool=false)
        !clone && return new(ptr)
        g = new(gclone(ptr))
        finalizer(g, gunclone)
        g
    end
end

function Gen(size::Integer)
    #    p = convert(GEN, pari_malloc(size*BYTES_IN_LONG))
    g = Gen(newblock(size))
    #    finalizer(g, g->pari_free(g.ptr))
    finalizer(g, g->killblock(g.ptr))
    g
end

Base.unsafe_convert(::Type{GEN}, x::Gen) = x.ptr

const GENGen = Union{GEN, Gen}



newblock(n::Integer) =
    ccall((:newblock, libpari), GEN, (Csize_t,), n)

GENtostr(x::GEN) = ccall((:GENtostr, libpari), Ptr{Cchar}, (GEN,), x)
gp_read_str(s) = ccall((:gp_read_str, libpari), GEN, (Ptr{Cchar},), s)




_libp(op::Symbol) = (op, :libpari)

for op in [:gunclone, :killblock]
    @eval begin
        $op(x::GENGen) = ccall($(_libp(op)), Void, (GEN,), x)
    end
end

for op in [:absi, :absr, :mpabs,
           :negi, :negr, :mpneg,
           :sqri, :sqrr, :mpsqr,
           :absi_shallow, # TODO: a "Gen" version probably doesn't make sense...
           :gclone]
    @eval begin
        $op(x::GEN) = ccall($(_libp(op)), GEN, (GEN,), x)
        $op(x::Gen) = keep_avma(()-> Gen($op(x.ptr), true))
    end
end

for op in [:mpcmp, :cmpii, :cmpir, :cmpri, :cmprr,
           :equalii, :equalrr,
           :absi_cmp, :absi_equal, :absr_cmp]
    @eval begin
        $op(x::G, y::G) where {G <: GENGen} =
            ccall($(_libp(op)), Cint, (GEN, GEN), x, y)
    end
end

for _op in [:add, :sub, :mul, :div, :rem, :mod],
    op in [Symbol(:mp, _op);
            [Symbol(_op, suf) for suf in [:ii, :ir, :ri, :rr]]]
    op in [:mprem, :mpmod] && continue
    @eval begin
        $op(x::GEN, y::GEN) = ccall($(_libp(op)), GEN, (GEN, GEN), x, y)
        $op(x::Gen, y::Gen) = keep_avma(()-> Gen($op(x.ptr, y.ptr), true))
    end
    opz = Symbol(op, :z)
    opz == :mpdivz && continue
    @eval begin
        $opz(x::GEN, y::GEN, z::GEN) =
            ccall($(_libp(opz)), GEN, (GEN, GEN, GEN), x, y, z)
        $opz(x::Gen, y::Gen, z::Gen) =
            keep_avma(()-> Gen($opz(x.ptr, y.ptr, z.ptr), true))
    end
end


## convenience functions

typ(x::Gen) = typ(x.ptr)

getavma() = unsafe_load(avma)

function __init__()
    pari_init(100_000_000, 2^16)
    atexit(pari_close)


    const global avma = cglobal((:avma, libpari), Ptr{Int})
    const global gen_0 = cglobal((:gen_0, libpari), Ptr{Int})
    const global gen_1 = cglobal((:gen_1, libpari), Ptr{Int})
    const global gen_2 = cglobal((:gen_2, libpari), Ptr{Int})
    const global gen_m1 = cglobal((:gen_m1, libpari), Ptr{Int})
    const global gen_m2 = cglobal((:gen_m2, libpari), Ptr{Int})

    const global pari_sigint = cglobal((:cb_pari_sigint, libpari), Ptr{Void})

    const global Gen_0 = Gen(unsafe_load(gen_0))
    const global Gen_1 = Gen(unsafe_load(gen_1))
    const global Gen_2 = Gen(unsafe_load(gen_2))
    const global Gen_M1 = Gen(unsafe_load(gen_m1))
    const global Gen_M2 = Gen(unsafe_load(gen_m2))

end


function pari_print(io::IO, a::Gen)
    cstr = GENtostr(a.ptr)
    print(io, unsafe_string(cstr))
    pari_free(cstr)
end

Base.show(io::IO, ::MIME{Symbol("text/plain")}, g::Gen) = pari_print(io, g)

function keep_avma(thunk)
    av = unsafe_load(avma)
    res = thunk()
    unsafe_store!(avma, av)
    res
end

macro pari_str(s)
    keep_avma(()-> Gen(gp_read_str(s), true))
end


mpodd(x::GENGen) = ccall((:mpodd, libpari), Cint, (GEN,), x)

type_name(t::Integer) = ccall((:type_name, libpari), Ptr{Cchar}, (long,), t)

type_name(x::Gen) = unsafe_string(type_name(typ(x.ptr)))

#### BigInt #########################
#=
function pari!(x::Ptr{Int}, a::fmpz, s::Int)
   siz = size(a)
   unsafe_store!(x, evaltyp(t_INT) | s, 1)
   unsafe_store!(x, evalsigne(sign(a)) | s, 2)
   z = BigInt(a)
   for i in 1:siz
      unsafe_store!(x, reinterpret(Int, unsafe_load(z.d, i)), i + 2)
   end
   return s
end
function pari(a::fmpz)
   s = gensize(a)
   g = pari_int(s)
   pari!(reinterpret(Ptr{Int}, g.d), a, s)
   return g
end
=#

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

Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::PariInt) = pari_print(io, x.i)



+(x::PariInt, y::PariInt) = PariInt(x.i + y.i)


end # module pari
