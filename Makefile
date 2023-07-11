CC = g++
# CC = clang
INCLUDES :=-I./pbis/include -I./
CFLAGS = -Wall -std=c++17 -O3
SOURCES = main.cpp

# will run pbis with sequential step, no execution parallel
# LDFLAGS := -fopenmp

# enables execution policy par (copy and other std algos)
LDFLAGS := -fopenmp -ltbb -D_USETBB

# will enable parallel mode for all algorithms
# LDFLAGS := -fopenmp -D_GLIBCXX_PARALLEL

unary:
	@echo " [BLD] Building binary unary. LDFLAGS: $(LDFLAGS)"
	@$(CC) $(CFLAGS) $(INCLUDES) $(SOURCES) $(LDFLAGS) -o prog
