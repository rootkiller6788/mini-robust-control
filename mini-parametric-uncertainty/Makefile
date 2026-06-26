# Makefile for mini-parametric-uncertainty
# Robust Control under Parametric Uncertainty

CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g -Wno-misleading-indentation -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable
INCLUDES = -I./include
LDFLAGS = -lm

SRC = src/pu_core.c src/pu_kharitonov.c src/pu_stability.c \
      src/pu_lmi.c src/pu_polytope.c src/pu_app.c
OBJ = $(SRC:.c=.o)
LIB = libpu.a

TEST_SRC = tests/test_core.c tests/test_kharitonov.c tests/test_stability.c \
           tests/test_lmi.c tests/test_polytope.c tests/test_app.c
TEST_BIN = $(TEST_SRC:.c=.exe)

EXAMPLE_SRC = examples/ex1_dc_motor.c examples/ex2_msd.c examples/ex3_kharitonov.c
EXAMPLE_BIN = $(EXAMPLE_SRC:.c=.exe)

.PHONY: all clean test examples lib

all: lib test examples

lib: $(LIB)

$(LIB): $(OBJ)
	ar rcs $@ $^

src/%.o: src/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

test: $(TEST_BIN)
	@echo "=== Running tests ==="
	@for t in $(TEST_BIN); do ./$$t || exit 1; done
	@echo "=== All tests passed ==="

tests/%.exe: tests/%.c $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L. -lpu -lm -o $@

examples: $(EXAMPLE_BIN)

examples/%.exe: examples/%.c $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L. -lpu -lm -o $@

clean:
	rm -f src/*.o $(LIB)
	rm -f tests/*.exe
	rm -f examples/*.exe
	rm -f _test.txt _bld.py _gen.py

check: test examples
	@echo "=== Running examples ==="
	@for e in $(EXAMPLE_BIN); do echo "--- $$e ---"; ./$$e; done

help:
	@echo "Targets: all, lib, test, examples, clean, check, help"
