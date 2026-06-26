# Makefile for mini-gap-metric-robustness
# Gap Metric and Nu-Gap Metric Robust Control Library

CC      ?= gcc
CFLAGS  ?= -std=c11 -Wall -Wextra -O2 -g
AR      ?= ar
ARFLAGS ?= rcs

INC_DIR  = include
SRC_DIR  = src
BLD_DIR  = build

INCLUDES = -I$(INC_DIR)
LDFLAGS  = -lm

SRCS = src/gap_core.c src/gap_system.c src/gap_coprime.c src/gap_metric.c src/gap_robustness.c src/gap_loopshape.c
OBJS = $(SRCS:src/%.c=$(BLD_DIR)/%.o)
LIB  = $(BLD_DIR)/libgap.a

.PHONY: all test clean

all: $(LIB)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(LIB): $(OBJS) | $(BLD_DIR)
	$(AR) $(ARFLAGS) $@ $^

$(BLD_DIR)/%.o: src/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

test: $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) tests/test_all.c -L$(BLD_DIR) -lgap $(LDFLAGS) -o $(BLD_DIR)/test_all
	./$(BLD_DIR)/test_all

examples: $(LIB)
	$(CC) $(CFLAGS) $(INCLUDES) examples/example_main.c -L$(BLD_DIR) -lgap $(LDFLAGS) -o $(BLD_DIR)/example_main
	./$(BLD_DIR)/example_main

clean:
	rm -rf $(BLD_DIR)
