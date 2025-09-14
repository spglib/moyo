.DEFAULT_GOAL := build

LOCAL_HOST := $(shell rustc -Vv | grep host | awk '{print $$2}')
TARGET_TRIPLE ?= $(LOCAL_HOST)

WORKSPACE_ROOT := $(shell cargo locate-project --workspace --message-format plain | xargs dirname)
VERSION := $(shell cargo metadata --no-deps --format-version 1 | jq -r '.packages[0].version')
DISTNAME := moyoc-$(VERSION)-$(TARGET_TRIPLE)
DISTDIR := $(WORKSPACE_ROOT)/moyoc/dist/$(DISTNAME)

# Infer dynamic library extension
ifeq ($(findstring windows,$(TARGET_TRIPLE)),windows)
    LIB_EXTENSION = dll
else ifeq ($(findstring apple,$(TARGET_TRIPLE)),apple)
    LIB_EXTENSION = dylib
else
    LIB_EXTENSION = so
endif

LIBNAME := libmoyoc.$(LIB_EXTENSION)
ifeq ($(TARGET_TRIPLE),$(LOCAL_HOST))
	LIBDIR := $(WORKSPACE_ROOT)/target/release
	CARGO_TARGET_ARG :=
else
	LIBDIR := $(WORKSPACE_ROOT)/target/$(TARGET_TRIPLE)/release
    CARGO_TARGET_ARG := --target $(TARGET_TRIPLE)
endif
BUILD_ARTIFACT := $(LIBDIR)/$(LIBNAME)

SRCS := $(wildcard src/*.rs src/**/*.rs)
TESTS := $(wildcard tests/*.c)
TEST_EXECUTABLES := $(TESTS:.c=.out)
TARGETS := $(LIBNAME) moyoc.h

libmoyoc.%: $(SRCS) Cargo.toml
	cargo build --release $(CARGO_TARGET_ARG)
	cp $(BUILD_ARTIFACT) .

moyoc.h: $(SRCS) cbindgen.toml
	cbindgen --config cbindgen.toml --crate moyoc --output moyoc.h

.PHONY: build
build: $(TARGETS)

.PHONY: dist
dist: build
	mkdir -p $(DISTDIR)/include $(DISTDIR)/lib
	cp moyoc.h $(DISTDIR)/include/
	cp $(LIBNAME) $(DISTDIR)/lib/

ifeq ($(LIB_EXTENSION),dll)
	cd dist && zip -r $(DISTNAME).zip $(DISTNAME)
else
	cd dist && tar czf $(DISTNAME).tar.gz $(DISTNAME)
endif
	@echo "Created distribution in dist/$(DISTNAME)"

tests/%.out: tests/%.c $(TARGETS)
	$(CC) $< -std=c11 -Wall -Wextra -I$(DISTDIR)/include -L$(DISTDIR)/lib -lmoyoc -fsanitize=address -fno-omit-frame-pointer -O2 -o $@

.PHONY: tests
tests: dist $(TEST_EXECUTABLES)
	@for test in $(TEST_EXECUTABLES); do \
		echo "Running $$test..."; \
		LD_LIBRARY_PATH=$(DISTDIR)/lib:$$LD_LIBRARY_PATH ./$$test; \
	done

.PHONY: clean
clean:
	rm -f libmoyoc.* moyoc.h $(TEST_EXECUTABLES)
	rm -rf dist
