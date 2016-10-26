
# Set the default VPATH
VPATH ?= $(shell pwd)

ARCH_MAKE_DEFAULT=$(VPATH)/arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)

.PHONY: default
default: all

.PHONY: all
all: lib
	@echo Done with everything

.PHONY: lib
lib:
ifdef LIBFDICT
	@echo "Using pre-built LIBFDICT: $(LIBFDICT)"
else
	$(MAKE) -C fdict "VPATH=$(VPATH)/fdict" \
		"ARCH_MAKE=$(ARCH_MAKE)" lib
endif
	$(MAKE) -C src "VPATH=$(VPATH)/src" lib

.PHONY: test
test: lib
	$(MAKE) -C test "VPATH=$(VPATH)/test" all

.PHONY: clean
clean:
ifndef LIBFDICT
	-if [ -d fdict ]; then $(MAKE) -C fdict "VPATH=$(VPATH)/fdict" \
		"ARCH_MAKE=$(ARCH_MAKE)" clean ; fi
endif
	$(MAKE) -C src clean
	$(MAKE) -C test clean
