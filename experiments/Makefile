# Compiler
CC := g++ -I./
NVCC := nvcc

# Files
COMMON_DIR := ./common/
COMMON_FILES := data_gen.cpp  utils.cpp sequential_impls.cpp 
COMMON_SRCS :=  $(addprefix $(COMMON_DIR), $(COMMON_FILES))

PARALLEL_FILES := parallel_impls.cpp experiments.cpp
PBIS_SRCS := Block-InsertionSort/pbis.cpp Block-InsertionSort/experimental_pbis.cpp
EXTERNAL_SRCS :=  ./external/mergesort/mergesort.cpp  
PARALLEL_SRCS := $(addprefix $(COMMON_DIR), $(PARALLEL_FILES)) $(PBIS_SRCS) $(EXTERNAL_SRCS)

MAIN_SRC := main.cpp

PARALLEL_OBJS := $(PARALLEL_SRCS:.cpp=.o)
COMMON_OBJS := $(COMMON_SRCS:.cpp=.o)
MAIN_OBJ := $(MAIN_SRC:.cpp=.o)
BIN := prog

# flags
INCLUDES :=-I./external/ips4o/include 
INCLUDES +=-I./external/ips2ra/include

# REG_SRC:= ./external/parallel-inplace-radixsort/regions_impl.cpp
# REG_INC :=-I./external/parallel-inplace-radixsort/RadixSort
# REG_INC +=-I./external/parallel-inplace-radixsort/common
# REG_OBJ := $(REG_SRC:.cpp=.o)

LDFLAGS :=-ltbb -fopenmp -latomic -lpthread 
# CFLAGS := -O3 -std=c++17 $(INCLUDES) -c
CFLAGS := -O3 -std=c++17 -Wall $(INCLUDES) -c
PARALLEL_FLAGS := -fopenmp -D_REENTRANT -DOPENMP

ALL_OBJS := $(patsubst %.cpp,%.o,$(COMMON_SRCS)) \
			$(patsubst %.cpp,%.o,$(SEQUENTIAL_SRCS)) \
			$(patsubst %.cpp,%.o,$(EXTERNAL_SRCS)) \
			$(patsubst %.cpp,%.o,$(PBIS_SRCS)) \
			$(patsubst %.cpp,%.o,$(PARALLEL_SRCS)) \
			$(patsubst %.cpp,%.o,$(MAIN_SRC)) \


# PBIS_SRCS, requires -fopenmp 
$(patsubst %.cpp,%.o,$(PARALLEL_SRCS)): %.o: %.cpp
	@echo "file: $@"
	$(CC) $(CFLAGS) $(PARALLEL_FLAGS) $(LDFLAGS) $< -o $@

$(patsubst %.cpp,%.o,$(COMMON_SRCS)): %.o: %.cpp
	@echo "file: $@"
	$(CC) $(CFLAGS) $< -o $@


# Build rule for main object file
$(MAIN_OBJ): $(MAIN_SRC)
	$(CC) $(CFLAGS) $< -o $@

# Rule to build the executable
prog: $(COMMON_OBJS) $(PARALLEL_OBJS) $(MAIN_OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)


# Make clean target phony
.PHONY: clean
# Clean rule
clean:
	rm -f $(BIN) $(ALL_OBJS)

# Remove clean target from default targets
.DEFAULT_GOAL := prog
