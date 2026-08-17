CC      ?= cc
PYTHON  ?= python3
BUILD   ?= build

WARN     := -Wall -Wextra -Wno-unused-parameter
OPT      := -O3
ifdef NATIVE
OPT      += -march=native
endif

ifneq ($(wildcard /usr/include/CL/cl.h),)
CL_INC   :=
else
CL_INC   := -Ithird_party/OpenCL-Headers
endif

UNAME_S  := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
CL_LIB   := -framework OpenCL
else
CL_LIB   := $(firstword $(wildcard \
              /usr/lib/$(shell uname -m)-linux-gnu/libOpenCL.so \
              /usr/lib/$(shell uname -m)-linux-gnu/libOpenCL.so.1 \
              /usr/lib64/libOpenCL.so /usr/lib64/libOpenCL.so.1 \
              /usr/lib/libOpenCL.so /usr/lib/libOpenCL.so.1))
endif

CFLAGS   := $(OPT) $(WARN) -I$(BUILD) $(CL_INC) $(EXTRA_CFLAGS)
LDLIBS   := $(CL_LIB) -lm

MONOCYPHER := third_party/monocypher/src/monocypher.c
CORE       := $(BUILD)/monocypher_x25519.inc
KERNELS    := $(BUILD)/mv_kernel_src.h
KERNEL_SRC := kernels/mv_prelude.h kernels/mv_edwards.h kernels/mv_grind.cl

.PHONY: all test format clean
all: mvgrind

$(BUILD):
	@mkdir -p $(BUILD)

$(MONOCYPHER):
	@echo "Monocypher submodule is missing. Run:"
	@echo "    git submodule update --init --recursive"
	@false

$(CORE): tools/extract_monocypher.py $(MONOCYPHER) | $(BUILD)
	$(PYTHON) tools/extract_monocypher.py $(MONOCYPHER) $@

$(KERNELS): tools/embed.py $(KERNEL_SRC) $(CORE) | $(BUILD)
	$(PYTHON) tools/embed.py $@ \
	  MV_SRC_PRELUDE=kernels/mv_prelude.h \
	  MV_SRC_X25519=$(CORE) \
	  MV_SRC_EDWARDS=kernels/mv_edwards.h \
	  MV_SRC_GRIND=kernels/mv_grind.cl

mvgrind: mvgrind.c $(KERNEL_SRC) $(CORE) $(KERNELS)
	@test -n "$(CL_LIB)" || { echo "no libOpenCL found: install your GPU vendor's OpenCL driver" >&2; exit 1; }
	$(CC) $(CFLAGS) -o $@ mvgrind.c $(LDLIBS)
	@echo "built: $@"

test: mvgrind
	./mvgrind --selftest

format:
	clang-format -i mvgrind.c $(KERNEL_SRC)

clean:
	rm -rf $(BUILD) mvgrind
