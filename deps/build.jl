push!(Libdl.DL_LOAD_PATH, "/usr/local/lib")

libpari = Libdl.dlopen(libpari)
libparimul = Libdl.dlsym(libpari, "__gmpn_mul")
if Integer(libparimul) == 0

end
