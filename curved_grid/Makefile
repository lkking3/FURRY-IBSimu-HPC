# Builds the curved-grid binary from ../sim/multi_grid_2d_curved.cpp
# Run `make` here or in ../sim/ (both produce the same binary).

SIM_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))../sim

.PHONY: all clean

all:
	$(MAKE) -C $(SIM_DIR) multi_grid_2d_curved

clean:
	$(MAKE) -C $(SIM_DIR) clean
