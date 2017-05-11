module Utils

export getavma, setavma!, keep_avma,
    at0, at, at1, at0!, at!, at1!, ul_,
    CulongMax, ClongMax

using Base.GMP: CulongMax, ClongMax

const GEN = Ptr{Clong}

at0(x) = unsafe_load(x)
at(x) = unsafe_load(x)
at(x, i) = unsafe_load(x, i+1)
at(x, i, j) = unsafe_load(unsafe_load(x, i+1), j+1)
at(x, i, j, k) = unsafe_load(unsafe_load(unsafe_load(x, i+1), j+1), k+1)
at1(x) = unsafe_load(x, 2)
at0!(x, a) = unsafe_store!(x, a)
at!(x, a) = unsafe_store!(x, a)
at1!(x, a) = unsafe_store!(x, a, 2)
ul_(x) = Base.unsafe_convert(Ptr{ulong}, x)


end # module Utils
