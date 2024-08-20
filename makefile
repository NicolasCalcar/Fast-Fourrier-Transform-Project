CC=gcc
CFLAGS=-Wall -lm 
OBJ=main.c

main: $(OBJ)
	$(CC) $(OBJ) -o main $(CFLAGS)

clean:
	rm main
