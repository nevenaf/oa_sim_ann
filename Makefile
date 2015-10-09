CC = gcc
CFLAGS = -Wall -std=c99 -fopenmp
DEPS = oa_search_functions.h
OBJ = run_oa_exchange.o oa_search_functions.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< -lm

run_oa_exchange: $(OBJ)
	gcc $(CFLAGS) -o $@ $^ -lm
