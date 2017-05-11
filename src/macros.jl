const LONG_IS_64_BIT = sizeof(Clong)*8 == 64 # defined from paricfg.h

# * from: parigen.h
if is_windows() && sizeof(Ptr{Void}) == 64 # TODO: _WIN64: check at runtime?
    const pari_ulong = Culonglong
    const long = Clonglong
else
    const pari_ulong = Culong
    const long = Clong
end
const ulong = pari_ulong

const GEN = Ptr{long}
const IUL = one(Culong) # == 1UL in C

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
const ULONG_MAX = ~(zero(Culong))

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


# * from: paricast.h

const _GEN = Ptr{GEN}
const __GEN = Ptr{Ptr{GEN}}
const ___GEN = Ptr{Ptr{Ptr{GEN}}}
const ____GEN = Ptr{Ptr{Ptr{Ptr{GEN}}}}
const _____GEN = Ptr{Ptr{Ptr{Ptr{Ptr{GEN}}}}}

mael2(m, i, j) = at(_GEN(m), i, j)
mael3(m, i, j, k) = at(__GEN(m), i, j, k)
mael4(m, i, j, k, l) = at(at(___GEN(m), i, j, k), l)
mael5(m, i, j, k, l, n) = at(at(____GEN(m), i, j, k), l, n)
const mael = mael2

gmael1(m, i) = at(_GEN(m), i)
gmael2(m, i, j) = at(__GEN(m), i, j)
gmael3(m, i, j, k) = at(___GEN(m), i, j, k)
gmael4(m, i, j, k, l) = at(at(____GEN(m), i, j, k), l)
gmael5(m, i, j, k, l, n) = at(at(_____GEN(m), i, j, k), l, n)
gmael(m, i, j) = gmael2(m, i, j)
gel(m, i) = gmael1(m, i)

gcoeff(a::GEN, i, j) = at(__GEN(a), j, i)
coeff(a::GEN, i, j) = at(_GEN(a), j, i)

GSTR(x) = Ptr{Cchar}(GEN(x) + sizeof(Clong))
