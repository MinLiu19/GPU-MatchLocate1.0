NVCC = nvcc -arch=sm_35  -Xcompiler "-fopenmp -fPIC" -Xlinker -lgomp 
CC = gcc -Os -mcmodel=medium -fopenmp
LIBS= -L/opt/sac/lib -lsac -lsacio -lm

BIN = ../bin

all: GPU_MatchLocate SelectFinal sac_e SHIFT lsac clean

GPU_MatchLocate: GPU_MatchLocate.cu sacio.o
	${NVCC} -o $(BIN)/$@ $^ $(LIBS)

SelectFinal: SelectFinal.o sacio.o
	${CC}   -o $(BIN)/$@ $^ $(LIBS)

sac_e: sac_e.o sacio.o
	${CC} 	-o $(BIN)/$@ $^ $(LIBS)

SHIFT: SHIFT.o sacio.o
	${CC}   -o $(BIN)/$@ $^ $(LIBS)

lsac: lsac.o sacio.o
	${CC}   -o $(BIN)/$@ $^ $(LIBS)

clean:
	-rm *.o
