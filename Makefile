CXXFLAGS = -O3 -Wall

.PHONY: all
all: blue-noise-dither blue-noise-ppm golden-ratio-rng blue-noise-process bayer-matrix

.PHONY: clean
clean:
	rm -f blue-noise-dither blue-noise-ppm golden-ratio-rng blue-noise-process bayer-matrix
