CXX = gcc
EXEC = main
CXXFLAGS = $(flag) -std=c++14 -Wall -pedantic
SRC = $(wildcard *.c)
OBJS = $(SRC:.c=.o)

all: $(EXEC) clear run

$(EXEC): $(OBJS)
	$(CXX) $^ -o $@

error: $(EXEC) run clear
error: flag=-D ERROR

clean:
	rm -rf *.o $(EXEC)

clear:
	rm -rf *.o

run:
	./$(EXEC)

check: 
	valgrind ./$(EXEC)

.PHONY: all clean clear run check