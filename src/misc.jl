using PARI.Lib


noexport = Dict(
PARI.Lib => Set(split("Lib eval at0 at at1 at0! at! at1! ul_"))
)

function make_exports(mod)
    exported = filter(s -> s[1] != '#' && !(s in noexport[mod]),
                      map(String, names(mod, true)))
    join(exported, ", ")
end
