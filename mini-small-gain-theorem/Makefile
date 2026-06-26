# Makefile for mini-small-gain-theorem
# Small Gain Theorem Module — Robust Control Theory

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g
LDFLAGS  = -lm
INC      = -I include

SRC_DIR  = src
TEST_DIR = tests
EX_DIR   = examples
BUILD_DIR = build

SRC = $(SRC_DIR)/sgt_core.c \
      $(SRC_DIR)/sgt_hinf.c \
      $(SRC_DIR)/sgt_stability.c \
      $(SRC_DIR)/sgt_uncertainty.c

OBJ = $(SRC:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

TEST_BIN = $(BUILD_DIR)/test_sgt
EX1_BIN  = $(BUILD_DIR)/example1_small_gain
EX2_BIN  = $(BUILD_DIR)/example2_dc_motor_robust
EX3_BIN  = $(BUILD_DIR)/example3_uncertainty

.PHONY: all test examples clean

all: $(OBJ) examples

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(TEST_BIN): $(OBJ) $(TEST_DIR)/test_sgt.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INC) $(TEST_DIR)/test_sgt.c $(OBJ) -o $@ $(LDFLAGS)

$(EX1_BIN): $(OBJ) $(EX_DIR)/example1_small_gain.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INC) $(EX_DIR)/example1_small_gain.c $(OBJ) -o $@ $(LDFLAGS)

$(EX2_BIN): $(OBJ) $(EX_DIR)/example2_dc_motor_robust.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INC) $(EX_DIR)/example2_dc_motor_robust.c $(OBJ) -o $@ $(LDFLAGS)

$(EX3_BIN): $(OBJ) $(EX_DIR)/example3_uncertainty.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INC) $(EX_DIR)/example3_uncertainty.c $(OBJ) -o $@ $(LDFLAGS)

test: $(TEST_BIN)
	@echo "=== Running Tests ==="
	$(TEST_BIN)

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)
	@echo "=== Examples built ==="

clean:
	rm -rf $(BUILD_DIR)
