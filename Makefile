CC ?= cc
CFLAGS ?= -O2 -Wall -Iinclude -std=c11
AR ?= ar
RANLIB ?= ranlib
PREFIX ?= /usr/local

SRC_DIR := src
INCLUDE_DIR := include
OBJ_DIR := build/obj
LIB_DIR := build/lib
BIN_DIR := build/bin
LIBNAME := dfx
LIB := $(LIB_DIR)/lib$(LIBNAME).a

SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJS := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRCS))

EXAMPLES_DIR := examples

.PHONY: all lib examples clean install

all: lib examples

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB): $(OBJS) | $(LIB_DIR)
	$(AR) rcs $@ $(OBJS)
	$(RANLIB) $@

lib: $(LIB)

examples: $(LIB) | $(BIN_DIR)
	$(MAKE) -C $(EXAMPLES_DIR) BUILD_BIN_DIR=$(abspath $(BIN_DIR)) LIB_DIR=$(abspath $(LIB_DIR)) CC="$(CC)" CFLAGS="$(CFLAGS) -I../include"

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

install: lib
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	mkdir -p $(DESTDIR)$(PREFIX)/include
	cp $(LIB) $(DESTDIR)$(PREFIX)/lib/
	cp $(INCLUDE_DIR)/*.h $(DESTDIR)$(PREFIX)/include/

clean:
	# Remove only build artifacts created by this Makefile, not other project files
	rm -rf $(OBJ_DIR) $(LIB_DIR) $(BIN_DIR)
CC ?= cc
CFLAGS ?= -O2 -Wall -Iinclude -std=c11
AR ?= ar
