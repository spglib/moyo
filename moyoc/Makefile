.DEFAULT_GOAL := build

SRCS = $(wildcard src/*.rs src/**/*.rs)
TESTS = $(wildcard tests/*.c)
TEST_EXECUTABLES = $(TESTS:.c=.out)

# Determine the dynamic library extension based on the platform
ifeq ($(OS),Windows_NT)
    LIB_EXTENSION = dll
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        LIB_EXTENSION = so
    endif
    ifeq ($(UNAME_S),Darwin)
        LIB_EXTENSION = dylib
    endif
endif
TARGETS = libmoyoc.$(LIB_EXTENSION) moyoc.h

libmoyoc.%: $(SRCS) Cargo.toml
	cargo build --release
	cp ../target/release/libmoyoc.$(LIB_EXTENSION) .

moyoc.h: $(SRCS) cbindgen.toml
	cbindgen --config cbindgen.toml --crate moyoc --output moyoc.h

.PHONY: build
build: $(TARGETS)

tests/%.out: tests/%.c $(TARGETS)
	$(CC) $< -std=c11 -Wall -Wextra -I. -L../target/release -lmoyoc -fsanitize=address -fno-omit-frame-pointer -O2 -g -o $@

.PHONY: test
test: build $(TEST_EXECUTABLES)
	@for test in $(TEST_EXECUTABLES); do \
		echo "Running $$test..."; \
		./$$test; \
	done

.PHONY: clean
clean:
	rm -f libmoyoc.* moyoc.h $(TEST_EXECUTABLES)
