module GP

export @gp_str, factor, matsize

import Base: length, setindex!, getindex

using PARI.Utils
using PARI.Lib

import PARI.Lib: typ, lg

typ(x::Gen) = typ(x.ptr)
lg(x::Gen) = lg(x.ptr)

function __init__()
#    pari_init_opts(100_000_000, 2^16, 4)
    pari_init(100_000_000, 2^16)
    atexit(pari_close)
end

macro gp_str(s)
    keep_avma(()-> Gen(gp_read_str(s), true))
end

factor(x::Gen, y::Gen = Gen_nil)::Gen = Lib.gp_factor0(x, y)

length(x::Gen) = Lib.glength(x)

getindex(x::Gen, i::Integer) = Gen(_getindex(x.ptr, Val{typ(x)}(), i))
getindex(x::Gen, i::Integer, j::Integer) =
    Gen(_getindex(x.ptr, Val{typ(x)}(), i, j))

_getindex(x::GEN, ::Union{Val{t_VEC}, Val{t_COL}, Val{t_MAT}}, i) = gel(x, i)
_getindex(x::GEN, ::Val{t_MAT}, i, j) = gcoeff(x, i, j)

end # module GP
