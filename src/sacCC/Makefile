CC = gcc
#Please change the direction to your SAC LIB.
CFLAGS = -L/usr/local/sac/lib -lsac -lsacio -lm  -fPIC -fno-pie -no-pie
all:ccsacc clean

BIN = ../../bin

ccsacc:ccsacc.o
	$(CC) -o $(BIN)/ccsacc ccsacc.o $(CFLAGS)

clean:
	rm -f *.o
