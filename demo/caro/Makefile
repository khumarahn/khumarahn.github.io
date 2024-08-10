OO = -I/usr/include/eigen3 -lembind --closure=1 -Wall -std=c++20 -sALLOW_MEMORY_GROWTH=1
DD = -I./eigen -Wall -std=c++20 -lgmp -lmpfr -fopenmp # -g -fsanitize=address -fsanitize=address -fsanitize=leak -fsanitize=undefined -ftrapv

.PHONY lsv::
	rm -f lsv-cpp.js lsv-cpp.wasm
	em++  $(OO) -O2 lsv.cpp -o lsv-cpp.js

.PHONY debug::
	g++  $(DD) -O3 -march=native lsv-debug.cpp -o lsv-cpp.x
