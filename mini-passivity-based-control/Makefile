# Makefile for mini-passivity-based-control
CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -O2 -g
INCLUDES = -Iinclude
LDFLAGS  = -lm

SRCDIR   = src
OBJDIR   = build
TESTDIR  = tests
EXDIR    = examples

SRCS = src/pbc_core.c src/pbc_passivity.c src/pbc_energy_shaping.c        src/pbc_port_hamiltonian.c src/pbc_ida.c src/pbc_el_systems.c        src/pbc_applications.c
OBJS = build/pbc_core.o build/pbc_passivity.o build/pbc_energy_shaping.o        build/pbc_port_hamiltonian.o build/pbc_ida.o build/pbc_el_systems.o        build/pbc_applications.o

.PHONY: all clean test

all: check

check:
	@mkdir -p $(OBJDIR)
	@echo "=== Compiling all sources ==="
	$(CC) $(CFLAGS) $(INCLUDES) -c src/pbc_core.c -o build/pbc_core.o $(LDFLAGS)
	$(CC) $(CFLAGS) $(INCLUDES) -c src/pbc_passivity.c -o build/pbc_passivity.o $(LDFLAGS)
	$(CC) $(CFLAGS) $(INCLUDES) -c src/pbc_energy_shaping.c -o build/pbc_energy_shaping.o $(LDFLAGS)
	$(CC) $(CFLAGS) $(INCLUDES) -c src/pbc_port_hamiltonian.c -o build/pbc_port_hamiltonian.o $(LDFLAGS)
	$(CC) $(CFLAGS) $(INCLUDES) -c src/pbc_ida.c -o build/pbc_ida.o $(LDFLAGS)
	$(CC) $(CFLAGS) $(INCLUDES) -c src/pbc_el_systems.c -o build/pbc_el_systems.o $(LDFLAGS)
	$(CC) $(CFLAGS) $(INCLUDES) -c src/pbc_applications.c -o build/pbc_applications.o $(LDFLAGS)
	@echo "=== All sources compiled ==="

test: check
	@mkdir -p $(OBJDIR)
	@for t in tests/*.c; do 		name=$$(basename $$t .c); 		echo "Building test: $$name"; 		$(CC) $(CFLAGS) $(INCLUDES) $$t $(OBJS) -o build/test_$$name $(LDFLAGS) || exit 1; 		echo "Running test_$$name..."; 		./build/test_$$name || exit 1; 	done
	@echo "=== All tests passed ==="

examples: check
	@mkdir -p $(OBJDIR)
	@for e in examples/*.c; do 		name=$$(basename $$e .c); 		echo "Building example: $$name"; 		$(CC) $(CFLAGS) $(INCLUDES) $$e $(OBJS) -o build/example_$$name $(LDFLAGS) || exit 1; 	done

clean:
	rm -rf $(OBJDIR)

lines:
	@echo "include/:"
	@wc -l include/*.h
	@echo "src/:"
	@wc -l src/*.c
	@echo -n "Total: "
	@cat include/*.h src/*.c | wc -l
