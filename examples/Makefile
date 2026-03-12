CC ?= cc
CFLAGS ?= -O2 -Wall -I../include
BUILD_BIN_DIR ?= ../build/bin
LIB_DIR ?= ../build/lib
LIB := $(LIB_DIR)/libdfx.a
LIBS ?= -lm

EXS := $(wildcard *.c)
BINS := $(patsubst %.c,$(BUILD_BIN_DIR)/%,$(EXS))

.PHONY: all clean

all: $(BINS)

$(BUILD_BIN_DIR)/%: %.c | $(BUILD_BIN_DIR)
	$(CC) $(CFLAGS) $< $(LIB) $(LIBS) -o $@

$(BUILD_BIN_DIR):
	mkdir -p $(BUILD_BIN_DIR)

clean:
	-rm -f $(BUILD_BIN_DIR)/*
