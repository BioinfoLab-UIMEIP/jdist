# jdist Makefile (macOS/Linux)
# - Builds src/jdist.cpp into ./jdist
# - Requires OpenCL headers + library and OpenMP support

CXX ?= g++
TARGET := jdist
SRC := src/jdist.cpp

# Common flags
CXXFLAGS ?= -O3 -std=c++17
WARNFLAGS := -Wall -Wextra -Wno-deprecated-declarations
OMPFLAG := -fopenmp

# Detect OS for OpenCL link flags
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  # macOS: OpenCL as a framework
  OPENCL_LDFLAGS := -framework OpenCL
else
  # Linux: OpenCL as a library
  OPENCL_LDFLAGS := -lOpenCL
endif

LDFLAGS := $(OPENCL_LDFLAGS)

.PHONY: all clean help run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(WARNFLAGS) $(OMPFLAG) -o $@ $< $(LDFLAGS)

run: $(TARGET)
	@echo "Running jdist on examples/example_matrix.tsv ..."
	./$(TARGET) examples/example_matrix.tsv examples/example_dist.tsv 4
	@echo "Wrote: examples/example_dist.tsv"

clean:
	rm -f $(TARGET) *.o

help:
	@echo "Targets:"
	@echo "  make        Build ./jdist"
	@echo "  make run    Run on examples/example_matrix.tsv"
	@echo "  make clean  Remove built binary"
