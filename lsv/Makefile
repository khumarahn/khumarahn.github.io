# things we need to build the code
LL = -I./eigen -Wall -std=c++23 -lgmp -lmpfr -lmpfi -fopenmp
# debug and safety checks (although these make the program slooooooooow...........)
#DD = -g -fsanitize=address,leak,undefined
# this is fassssster
#DD = -g -fsanitize=undefined
# this is for profiling
DD = -g -fno-omit-frame-pointer

.PHONY rm.lsv-cpp.x::
	rm -f lsv-cpp.x

debug: lsv-cpp.x
lsv-cpp.x: *.h *.cpp Makefile
	clang++ $(LL) $(DD) -O3 -march=native lsv-debug.cpp -o lsv-cpp.x

loop: rm.lsv-cpp.x debug
	time ./lsv-cpp.x && true

.PHONY test::
	clang++ $(LL) $(DD) -O3 -march=native lsv-test.cpp -o lsv-test.x

.PHONY web::
	em++ -O2 -std=c++23 -flto -msimd128 \
		-Wall -Wno-deprecated-declarations \
		-lembind \
		-I../../../eigen \
		-I../../../boost_1_90_0 \
		-I../../../gmp-6.3.0/wasm-build/include \
		-I../../../mpfr-4.2.2/wasm-build/include \
		-I../../../mpfi-1.5.4/wasm-build/include \
		lsv.cpp \
		../../../mpfi-1.5.4/wasm-build/lib/libmpfi.a \
		../../../mpfr-4.2.2/wasm-build/lib/libmpfr.a \
		../../../gmp-6.3.0/wasm-build/lib/libgmp.a \
		-o lsv-cpp.js \
		-s ALLOW_MEMORY_GROWTH=1 -s STACK_SIZE=67108864
