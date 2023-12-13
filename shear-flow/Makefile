HH = -I../hardballs/src

.PHONY 20::
	rm -f 20.js 20.wasm
	em++ -lembind -s ALLOW_MEMORY_GROWTH=1 -Wall --closure=1 -O3 $(HH) 20.cpp -o 20.js
