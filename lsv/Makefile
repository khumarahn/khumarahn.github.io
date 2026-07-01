# things we need to build the code
LL = -I../external/eigen -Wall -std=c++23 -lgmp -lmpfr -lmpfi -fopenmp
# debug and safety checks (although these make the program slooooooooow...........)
#DD = -g -fsanitize=address,leak,undefined
# this is fassssster
DD = -g -fsanitize=undefined
# this is for profiling
#DD = -g -fno-omit-frame-pointer

.PHONY test::
	clang++ $(LL) $(DD) -O3 -march=native lsv-test.cpp -o lsv-test.x

.PHONY web::
	em++ -O2 -std=c++23 -flto -msimd128 \
		-Wall -Wno-deprecated-declarations \
		-lembind \
		-I../external/eigen \
		-I../external/gmp/wasm-build/include \
		-I../external/mpfr/wasm-build/include \
		-I../external/mpfi/wasm-build/include \
		lsv.cpp \
		../external/mpfi/wasm-build/lib/libmpfi.a \
		../external/mpfr/wasm-build/lib/libmpfr.a \
		../external/gmp/wasm-build/lib/libgmp.a \
		-o lsv-cpp.js \
		-s USE_BOOST_HEADERS=1 \
		-s ALLOW_MEMORY_GROWTH=1 -s STACK_SIZE=67108864
