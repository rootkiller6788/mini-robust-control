# Makefile for mini-mu-synthesis
# u-Synthesis: Structured Singular Value and DK-Iteration
# Reference: Packard & Doyle (1993), Zhou, Doyle & Glover (1996)

CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -O2
LDFLAGS  = -lm
INCLUDES = -Iinclude

SRCDIR   = src
TESTDIR  = tests
EXAMDIR  = examples
BLDDIR   = build

SRCS     = $(wildcard $(SRCDIR)/*.c)
OBJS     = $(patsubst $(SRCDIR)/%.c,$(BLDDIR)/%.o,$(SRCS))
TESTS    = $(wildcard $(TESTDIR)/*.c)
TESTBINS = $(patsubst $(TESTDIR)/%.c,$(BLDDIR)/test_%,$(TESTS))
EXAMPLES = $(wildcard $(EXAMDIR)/*.c)
EXAMBINS = $(patsubst $(EXAMDIR)/%.c,$(BLDDIR)/example_%,$(EXAMPLES))

.PHONY: all test examples clean

all: $(BLDDIR)/libmu.a $(TESTBINS) $(EXAMBINS)

$(BLDDIR):
	mkdir -p $(BLDDIR)

$(BLDDIR)/%.o: $(SRCDIR)/%.c | $(BLDDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(BLDDIR)/libmu.a: $(OBJS)
	ar rcs $@ $(OBJS)

$(BLDDIR)/test_%: $(TESTDIR)/%.c $(BLDDIR)/libmu.a | $(BLDDIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(BLDDIR) -lmu $(LDFLAGS) -o $@

$(BLDDIR)/example_%: $(EXAMDIR)/%.c $(BLDDIR)/libmu.a | $(BLDDIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(BLDDIR) -lmu $(LDFLAGS) -o $@

test: $(TESTBINS)
	@for t in $(TESTBINS); do echo "=== Running $$t ==="; $$t || exit 1; done
	@echo "All tests passed."

examples: $(EXAMBINS)
	@for e in $(EXAMBINS); do echo "=== Example $$e ==="; $$e || exit 1; done

clean:
	rm -rf $(BLDDIR)
