CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17
PKG   := ibsimu-1.0.6dev

all: two_grid_2d

two_grid_2d: two_grid_2d.cpp
	$(CXX) $(CXXFLAGS) two_grid_2d.cpp \
	  $(shell pkg-config --cflags --libs $(PKG)) \
	  -o $@

clean:
	rm -f two_grid_2d two_grid_geom.png

