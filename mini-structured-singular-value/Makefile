# Makefile for mini-structured-singular-value
# Structured Singular Value (mu) Analysis and Synthesis Library
#
# Build: make
# Test:  make test
# Clean: make clean

CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c11 -Iinclude
LDFLAGS = -lm

# Source files
SRC = src/ssv_core.c src/ssv_bounds.c src/ssv_uncertainty.c src/ssv_lft.c src/ssv_synthesis.c
OBJ = $(SRC:.c=.o)

# Test files
TEST_SRC = tests/test_ssv.c
TEST_BIN = test_runner

# Example files
EXAMPLES = examples/ex1_mu_analysis examples/ex2_uncertainty examples/ex3_dk_synthesis

.PHONY: all clean test examples

all: libssv.a

libssv.a: $(OBJ)
	ar rcs $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) -o $@ $(TEST_SRC) $(OBJ) $(LDFLAGS)

examples: $(EXAMPLES)

examples/ex1_mu_analysis: examples/ex1_mu_analysis.c $(OBJ)
	$(CC) $(CFLAGS) -o $@ $< $(OBJ) $(LDFLAGS)

examples/ex2_uncertainty: examples/ex2_uncertainty.c $(OBJ)
	$(CC) $(CFLAGS) -o $@ $< $(OBJ) $(LDFLAGS)

examples/ex3_dk_synthesis: examples/ex3_dk_synthesis.c $(OBJ)
	$(CC) $(CFLAGS) -o $@ $< $(OBJ) $(LDFLAGS)

clean:
	rm -f src/*.o $(TEST_BIN) libssv.a examples/ex1_mu_analysis examples/ex2_uncertainty examples/ex3_dk_synthesis
