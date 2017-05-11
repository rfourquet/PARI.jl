module Lib

using PARI.Utils
using Base: unsafe_convert

include("macros.jl")






# Warning: for native kernel
int_W(x, l) = x+(lgefint(x)-1-l)*BYTES_IN_LONG
int_W_lg(x, l, lx) = x+(lx-1-l)*BYTES_IN_LONG


# end macros

const libpari = :libpari

#push!(Libdl.DL_LOAD_PATH, "/usr/local/lib")
#Libdl.dlopen(libpari)

pari_init(size::Integer, maxprime::Integer) =
    ccall((:pari_init, libpari), Void, (Csize_t, Culong), size, maxprime)

pari_init_opts(size::Integer, maxprime::Integer, opts) =
    ccall((:pari_init, libpari), Void, (Csize_t, Culong, Culong),
          size, maxprime, opts)

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
        export $op
    end
end

for op in [:absi, :absr, :mpabs,
           :negi, :negr, :mpneg,
           :sqri, :sqrr, :mpsqr,
           :absi_shallow, # TODO: a "Gen" version probably doesn't make sense...
           :gclone,
           :Z_factor]
    @eval begin
        $op(x::GEN) = ccall($(_libp(op)), GEN, (GEN,), x)
        $op(x::Gen) = keep_avma(()-> Gen($op(x.ptr), true))
        export $op
    end
end

for op in [:glength]
    @eval begin
        $op(x::GEN) = ccall($(_libp(op)), Clong, (GEN,), x)
        $op(x::Gen) = $op(x.ptr)
        export $op
    end
end


for op in [:gp_factor0]
    @eval begin
        $op(x::GEN, y::GEN) = ccall($(_libp(op)), GEN, (GEN, GEN), x, y)
        $op(x::Gen, y::Gen) = keep_avma(()-> Gen($op(x.ptr, y.ptr), true))
        export $op
    end
end

for op in [:mpcmp, :cmpii, :cmpir, :cmpri, :cmprr,
           :equalii, :equalrr,
           :absi_cmp, :absi_equal, :absr_cmp]
    @eval begin
        $op(x::G, y::G) where {G <: GENGen} =
            ccall($(_libp(op)), Cint, (GEN, GEN), x, y)
        export $op
    end
end

for _op in [:add, :sub, :mul, :div, :rem, :mod],
    op in [Symbol(:mp, _op);
            [Symbol(_op, suf) for suf in [:ii, :ir, :ri, :rr]]]
    op in [:mprem, :mpmod] && continue
    @eval begin
        $op(x::GEN, y::GEN) = ccall($(_libp(op)), GEN, (GEN, GEN), x, y)
        $op(x::Gen, y::Gen) = keep_avma(()-> Gen($op(x.ptr, y.ptr), true))
        export $op
    end
    opz = Symbol(op, :z)
    opz == :mpdivz && continue
    @eval begin
        $opz(x::GEN, y::GEN, z::GEN) =
            ccall($(_libp(opz)), GEN, (GEN, GEN, GEN), x, y, z)
        $opz(x::Gen, y::Gen, z::Gen) =
            keep_avma(()-> Gen($opz(x.ptr, y.ptr, z.ptr), true))
        export $opz
    end
end

pari_sigint_handler() = error("User interrupt")
pari_handle_exception(x) = (println("ERROR: ", x); sleep(1);0 % Cint)

global _INITED = false
function __init__()
    global _INITED
    _INITED && return # avoid multiple initializations
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

    const global Gen_nil = Gen(GEN(0))

    const global pari_sigint = cglobal((:cb_pari_sigint, libpari), Ptr{Void})
    const global cb_pari_handle_exception = cglobal((:cb_pari_handle_exception, libpari), Ptr{Void})
   unsafe_store!(pari_sigint, cfunction(pari_sigint_handler, Void, ()), 1)
   unsafe_store!(cb_pari_handle_exception,
                cfunction(pari_handle_exception, Cint, (Clong,)), 1)
    _INITED = true
end


getavma() = unsafe_load(avma)
setavma!(av::GEN) = unsafe_store!(avma, av)

function keep_avma(thunk)
    av = getavma()
    res = thunk()
    setavma!(av)
    res
end


mpodd(x::GENGen) = ccall((:mpodd, libpari), Cint, (GEN,), x)

type_name(t::Integer) = ccall((:type_name, libpari), Ptr{Cchar}, (long,), t)

type_name(x::Gen) = unsafe_string(type_name(typ(x.ptr)))


export BIGDEFAULTPREC, BITS_IN_HALFULONG, BITS_IN_LONG, BYTES_IN_LONG,
    CLONEBIT, CLOSURE, COL, COMPLEX, DEFAULTPREC, ERROR, EXPOBITS,
    EXPOnumBITS, EXTRAPRECWORD, FFELT, FRAC, GEN, GENGen, GENtostr,
    Gen, Gen_0, Gen_1, Gen_2, Gen_M1, Gen_M2, HIGHBIT, HIGHEXPOBIT,
    HIGHMASK, HIGHVALPBIT, HIGHWORD, INT, INTMOD, IUL, LGBITS,
    LGnumBITS, LIST, LONG_IS_64_BIT, LONG_MAX, LOWDEFAULTPREC,
    LOWMASK, LOWWORD, MAT, MAXVARN, MEDDEFAULTPREC, NO_VARIABLE,
    PADIC, POL, POLMOD, PRECPBITS, PRECPSHIFT, QFI, QFR, QUAD, REAL,
    RFRAC, SER, SIGNBITS, SIGNSHIFT, SIGNnumBITS, SMALL_ULONG, STR,
    TWOPOTBITS_IN_LONG, TYPBITS, TYPSHIFT, TYPnumBITS, Types,
    ULONG_MAX, VALPBITS, VALPnumBITS, VARNBITS, VARNSHIFT,
    VARNnumBITS, VEC, VECSMALL, Z_factor, __init__, _evalexpo,
    _evallg, _evalprecp, _evalvalp, _libp, absi, absi_cmp, absi_equal,
    absi_shallow, absr, absr_cmp, addii, addiiz, addir, addirz, addri,
    addriz, addrr, addrrz, avma, cmpii, cmpir, cmpri, cmprr, divii,
    diviiz, divir, divirz, divri, divriz, divrr, divrrz, equalii,
    equalrr, evallg, evallgefint, evallgeflist, evalsigne,
    evaltyp, evalvarn, gclone, gen_0, gen_1, gen_2, gen_m1, gen_m2,
    getavma, gp_read_str, gunclone, int_W, int_W_lg, isclone,
    keep_avma, killblock, lg, lgefint, libpari, long, modii, modiiz,
    modir, modirz, modri, modriz, modrr, modrrz, mpabs, mpadd, mpaddz,
    mpcmp, mpdiv, mpmul, mpmulz, mpneg, mpodd, mpsqr, mpsub, mpsubz,
    mulii, muliiz, mulir, mulirz, mulri, mulriz, mulrr, mulrrz, negi,
    negr, newblock, pari_close, pari_free, pari_init, pari_malloc,
    pari_sigint, pari_ulong, remii, remiiz, remir, remirz, remri,
    remriz, remrr, remrrz, setavma!, setisclone, setlg, setlgefint,
    setsigne, settyp, signe, sqri, sqrr, subii, subiiz, subir, subirz,
    subri, subriz, subrr, subrrz, t_CLOSURE, t_COL, t_COMPLEX,
    t_ERROR, t_FFELT, t_FRAC, t_INT, t_INTMOD, t_LIST, t_MAT, t_PADIC,
    t_POL, t_POLMOD, t_QFI, t_QFR, t_QUAD, t_REAL, t_RFRAC, t_SER,
    t_STR, t_VEC, t_VECSMALL, typ, type_name, ulong, unsetisclone,
    pari_init_opts, Gen_nil, gel, gcoeff

end # module Lib
