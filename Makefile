.PHONY all::
	emcc -v >> emcc.txt 2>&1
	$(MAKE) -C demo
	$(MAKE) -C bb
	$(MAKE) -C shear-flow
