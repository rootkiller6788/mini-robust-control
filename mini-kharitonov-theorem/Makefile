# ==============================================================
# Makefile - mini-kharitonov-theorem
# Kharitonov Theorem for Robust Polynomial Stability
# ==============================================================

CC       = gcc
CFLAGS   = -Wall -Wextra -O2 -std=c99 -Iinclude
LDFLAGS  = -lm

SRCDIR   = src
INCDIR   = include
TESTDIR  = tests
EXAMDIR  = examples
BUILDDIR = build

CORE_SRC = $(SRCDIR)/kh_core.c $(SRCDIR)/kh_interval.c $(SRCDIR)/kh_hurwitz.c $(SRCDIR)/kh_construction.c $(SRCDIR)/kh_verification.c $(SRCDIR)/kh_application.c
CORE_OBJ = $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(CORE_SRC))

TEST_SRC = $(TESTDIR)/test_kharitonov.c
TEST_BIN = $(BUILDDIR)/test_kharitonov

EXAMPLE_BINS = $(BUILDDIR)/example_dc_motor $(BUILDDIR)/example_interval_robustness $(BUILDDIR)/example_aircraft_pitch

.PHONY: all test clean examples run-examples help

all: $(BUILDDIR) examples

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(CORE_OBJ) $(TEST_SRC) | $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $(TEST_SRC) $(CORE_OBJ) $(LDFLAGS)

test: $(TEST_BIN)
	$(TEST_BIN)

$(BUILDDIR)/example_dc_motor: $(CORE_OBJ) $(EXAMDIR)/example_dc_motor.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $(EXAMDIR)/example_dc_motor.c $(CORE_OBJ) $(LDFLAGS)

$(BUILDDIR)/example_interval_robustness: $(CORE_OBJ) $(EXAMDIR)/example_interval_robustness.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $(EXAMDIR)/example_interval_robustness.c $(CORE_OBJ) $(LDFLAGS)

$(BUILDDIR)/example_aircraft_pitch: $(CORE_OBJ) $(EXAMDIR)/example_aircraft_pitch.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $(EXAMDIR)/example_aircraft_pitch.c $(CORE_OBJ) $(LDFLAGS)

examples: $(EXAMPLE_BINS)

run-examples: examples
	$(BUILDDIR)/example_dc_motor
	$(BUILDDIR)/example_interval_robustness
	$(BUILDDIR)/example_aircraft_pitch

clean:
	rm -rf $(BUILDDIR)

help:
	@echo "Targets: all test examples clean run-examples help"
