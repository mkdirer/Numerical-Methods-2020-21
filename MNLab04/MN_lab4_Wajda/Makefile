CC=g++
CFLAGS= -O2 -lm -lm -lgsl -lgslcblas
DEPS=
OBJS=main.o
OUT=main

$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CFLAGS)
	
.c.o:
	$(CC) -c -o $@ $< $(CFLAGS) 
	
.PHONY: clean	
clean:
	rm -f *o
	rm -f *~
	rm -f $(OUT)
	
run:
	./main