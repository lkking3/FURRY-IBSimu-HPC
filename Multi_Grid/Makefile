# Binaries are built from ../sim/ — run `make` there first.
# This Makefile is kept for convenience: it delegates to sim/.

SIM_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))../sim

.PHONY: all clean

all:
	$(MAKE) -C $(SIM_DIR) multi_grid_2d

clean:
	$(MAKE) -C $(SIM_DIR) clean
