OBJS 	= ergasia3.o previousFunctions.o 
SOURCE	= ergasia3.c previousFunctions.c 
HEADER  = previousFunctions.h
OUT  	= ergasia3
CC	= gcc
CFLAGS   = -g -c 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

ergasia3.o: ergasia3.c
	$(CC) $(CFLAGS) ergasia3.c

previousFunctions.o: previousFunctions.c
	$(CC) $(CFLAGS) previousFunctions.c

clean:
	rm -f $(OBJS) $(OUT)

count:
	wc $(SOURCE) $(HEADER)