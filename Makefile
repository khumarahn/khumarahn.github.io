.PHONY all::
	emcc -v >> emcc.txt 2>&1
	$(MAKE) -C bb
