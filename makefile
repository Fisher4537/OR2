#all: mainTSP

mainTSP: mainTSP.o
	gcc -o mainTSP mainTSP.o

mainTSP.o: mainTSP.c
	gcc -c mainTSP.c

run:
	./mainTSP $(ARGS) >> log/out_`data`.txt

debug:
	@echo "compiling with debugging options..."
	gcc -g mainTSP.c -o mainTSP.x

clean:
	rm mainTSP.o mainTSP
