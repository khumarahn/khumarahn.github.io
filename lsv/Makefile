.DEFAULT_GOAL := web

WASM_OPT = -O3 -flto
GMP_LIB  = ../external/gmp/wasm-build/lib/libgmp.a
MPFR_LIB = ../external/mpfr/wasm-build/lib/libmpfr.a
MPFI_LIB = ../external/mpfi/wasm-build/lib/libmpfi.a

$(GMP_LIB):
	cd ../external/gmp && touch doc/version.texi && \
	(test -f configure || ./autogen.sh || autoreconf -fi) && \
	emconfigure ./configure --disable-assembly --host=none --enable-cxx --prefix=$$(pwd)/wasm-build CFLAGS="$(WASM_OPT)" CXXFLAGS="$(WASM_OPT)" && \
	emmake make -j$$(nproc) && emmake make install

$(MPFR_LIB): $(GMP_LIB)
	cd ../external/mpfr && \
	(test -f configure || ./autogen.sh || autoreconf -fi) && \
	emconfigure ./configure --host=none --with-gmp=$$(pwd)/../gmp/wasm-build --prefix=$$(pwd)/wasm-build CFLAGS="$(WASM_OPT)" CXXFLAGS="$(WASM_OPT)" && \
	emmake make -j$$(nproc) && emmake make install

$(MPFI_LIB): $(MPFR_LIB)
	cd ../external/mpfi && \
	(test -f configure || ./autogen.sh || autoreconf -fi) && \
	emconfigure ./configure --host=none --with-gmp=$$(pwd)/../gmp/wasm-build --with-mpfr=$$(pwd)/../mpfr/wasm-build --prefix=$$(pwd)/wasm-build CFLAGS="$(WASM_OPT)" CXXFLAGS="$(WASM_OPT)" && \
	emmake make -j$$(nproc) && emmake make install

# Note: web target now depends on the three library files and uses -O3
.PHONY web:: $(GMP_LIB) $(MPFR_LIB) $(MPFI_LIB)
	em++ -O3 -std=c++23 -flto -msimd128 \
		-Wall -Wno-deprecated-declarations \
		-lembind \
		-I../external/eigen \
		-I../external/gmp/wasm-build/include \
		-I../external/mpfr/wasm-build/include \
		-I../external/mpfi/wasm-build/include \
		lsv.cpp \
		$(MPFI_LIB) \
		$(MPFR_LIB) \
		$(GMP_LIB) \
		-o lsv-cpp.js \
		-s USE_BOOST_HEADERS=1 \
		-s ALLOW_MEMORY_GROWTH=1 -s STACK_SIZE=67108864

# LOCAL BUILD

LL = -I../external/eigen -Wall -std=c++23 -lgmp -lmpfr -lmpfi -fopenmp

#DD = -g -fsanitize=address,leak,undefined
DD = -g -fsanitize=undefined

.PHONY test::
	clang++ $(LL) $(DD) -O3 -march=native lsv-test.cpp -o lsv-test.x
