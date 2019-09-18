BIN=../../bin

all: sac_e

sac_e: %:%.o sacio.o
	$(LINK.c) -o ${BIN}/$@ $@.o sacio.o -lm	


clear:
	rm *.o
