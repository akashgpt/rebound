export OPENGL=0
export MPI=1
export CC=mpicc
include ../../src/Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -g -I /home/akashgpt/external_libraries/gsl/include/  -I ../../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) minor-planet_ring_parallel.c -L /home/akashgpt/external_libraries/gsl/lib/ -L.  -lrebound $(LIB) -lgsl -lgslcblas -lm -pthread -o rebound
	@echo ""
	@echo "REBOUND compiled successfully."

librebound: 
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../../src/
	@-rm -f librebound.so
	@ln -s ../../src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
