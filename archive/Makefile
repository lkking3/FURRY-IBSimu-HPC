CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17
PKG   := ibsimu-1.0.6dev

BINARIES := two_grid_2d multi_grid_2d

all: $(BINARIES)

two_grid_2d: two_grid_2d.cpp
	$(CXX) $(CXXFLAGS) two_grid_2d.cpp 	  $(shell pkg-config --cflags --libs $(PKG)) 	  -o $@

multi_grid_2d: multi_grid_2d.cpp
	$(CXX) $(CXXFLAGS) multi_grid_2d.cpp 	  $(shell pkg-config --cflags --libs $(PKG)) 	  -o $@

clean:
	rm -f $(BINARIES) two_grid_geom.png multi_grid_geom.png
