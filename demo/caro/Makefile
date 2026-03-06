# compilation for webassembly
OO = -I/usr/include/eigen3 -lembind --closure=1 -Wall -std=c++23 -sALLOW_MEMORY_GROWTH=1
# things we need to build the code
LL = -I./eigen -Wall -std=c++23 -lgmp -lmpfr -lmpfi -fopenmp
# debug and safety checks (although these make the program slooooooooow...........)
DD = -g -fsanitize=address,leak,undefined
# this is fassssster
#DD = -g -fsanitize=undefined

.PHONY lsv::
	rm -f lsv-cpp.js lsv-cpp.wasm
	em++  $(OO) -O2 lsv.cpp -o lsv-cpp.js

debug: lsv-cpp.x
lsv-cpp.x: *.h *.cpp Makefile
	clang++ $(LL) $(DD) -O3 -march=native lsv-debug.cpp -o lsv-cpp.x

.PHONY test::
	clang++ $(LL) $(DD) -O3 -march=native test.cpp -o test.x
