HH = -I/usr/include/eigen3

.PHONY em::
	rm -f em.js em.wasm
	em++ -lembind -Wall $(HH) -O2 em.cpp -o em.js
