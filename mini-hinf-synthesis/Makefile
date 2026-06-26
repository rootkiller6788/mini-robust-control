# mini-hinf-synthesis Makefile
# H-infinity Controller Synthesis Library
#
# Targets:
#   make          - build static library libhinf.a
#   make test     - build and run all tests
#   make examples - build all examples
#   make clean    - remove build artifacts

CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -O2 -g
AR       = ar
ARFLAGS  = rcs

INC_DIR  = include
SRC_DIR  = src
TST_DIR  = tests
EXA_DIR  = examples
BCH_DIR  = benches
BLD_DIR  = build

INCLUDES = -I$(INC_DIR)
LDFLAGS  = -lm

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(BLD_DIR)/%.o, $(SRCS))
LIB  = $(BLD_DIR)/libhinf.a

# Test files
TESTS_SRC = $(wildcard $(TST_DIR)/*.c)
TESTS_BIN = $(patsubst $(TST_DIR)/%.c, $(BLD_DIR)/%, $(TESTS_SRC))

# Example files
EXAMPLES_SRC = $(wildcard $(EXA_DIR)/*.c)
EXAMPLES_BIN = $(patsubst $(EXA_DIR)/%.c, $(BLD_DIR)/%, $(EXAMPLES_SRC))

# Bench files
BENCHES_SRC = $(wildcard $(BCH_DIR)/*.c)
BENCHES_BIN = $(patsubst $(BCH_DIR)/%.c, $(BLD_DIR)/%, $(BENCHES_SRC))

.PHONY: all test examples bench clean

all: $(LIB)

# Build library
$(LIB): $(OBJS)
	@mkdir -p $(BLD_DIR)
	$(AR) $(ARFLAGS) $@ $^

# Compile source files
$(BLD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Compile and run tests
test: $(LIB) $(TESTS_BIN)
	@echo "=== Running tests ==="
	@passed=0; failed=0; \
	for t in $(TESTS_BIN); do \
		if ./$$t; then passed=$$((passed+1)); \
		else failed=$$((failed+1)); fi; \
	done; \
	echo "=== $$passed passed, $$failed failed ==="

# Compile test binaries
$(BLD_DIR)/%: $(TST_DIR)/%.c $(LIB)
	@mkdir -p $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(BLD_DIR) -lhinf $(LDFLAGS) -o $@

# Build examples
examples: $(LIB) $(EXAMPLES_BIN)

$(BLD_DIR)/%: $(EXA_DIR)/%.c $(LIB)
	@mkdir -p $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(BLD_DIR) -lhinf $(LDFLAGS) -o $@

# Build benchmarks
bench: $(LIB) $(BENCHES_BIN)

$(BLD_DIR)/%: $(BCH_DIR)/%.c $(LIB)
	@mkdir -p $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(BLD_DIR) -lhinf $(LDFLAGS) -o $@

# Clean
clean:
	rm -rf $(BLD_DIR)
