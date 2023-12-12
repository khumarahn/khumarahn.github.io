HH = -I../hardballs/src

.PHONY 32::
	rm -f 32.js 32.wasm
	em++ -lembind -s ALLOW_MEMORY_GROWTH=1 -Wall --closure=1 -O3 $(HH) 32.cpp -o 32.js
