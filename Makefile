CC = g++
INCLUDES :=-I./pbis/include -I./
CFLAGS = -Wall -std=c++17 -O3  
SOURCES = main.cpp
LDFLAGS := -fopenmp -ltbb -D_USETBB
# LDFLAGS := -fopenmp

unary:
	@echo " [BLD] Building binary unary"
	@$(CC) $(CFLAGS) $(INCLUDES) $(SOURCES) $(LDFLAGS) -o prog
