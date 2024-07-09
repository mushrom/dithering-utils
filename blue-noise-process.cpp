#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdint>

#include <vector>
#include <array>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

struct image {
	uint8_t *data;
	int width;
	int height;
	int channels;
	size_t pixels;
	size_t length;
};

image load_image(const char *filename) {
	image ret;
	ret.data = stbi_load(filename, &ret.width, &ret.height, &ret.channels, 0);
	ret.pixels = ret.width * ret.height;
	ret.length = ret.pixels * ret.channels;

	printf("%s: %p channels=%d width=%d height=%d pixels=%lu length=%lu\n",
	       filename, ret.data, ret.channels, ret.width, ret.height, ret.pixels, ret.length);
	return ret;
}

int main(int argc, char *argv[]) {
	srand(time(NULL));
	if (argc < 4) {
		fprintf(stderr, "usage: ./process [input] [blue noise] [output]\n");
		return 1;
	}

	image in    = load_image(argv[1]);
	image noise = load_image(argv[2]);
	FILE *out = fopen(argv[3], "w");

	fprintf(out, "P4\n%d %d\n", in.width, in.height);

	unsigned bits = 0;
	int buffer = 0;

	for (int y = 0; y < in.height; y++) {
		for (int x = 0; x < in.width; x++) {
			size_t k = ((y*in.width) + x) * in.channels;
			size_t m = (((y % noise.height)*noise.width) + (x % noise.width)) * noise.channels;
			//float r = pow(in.data[k+0] / 255.0, 2.2);
			//float g = pow(in.data[k+1] / 255.0, 2.2);
			//float b = pow(in.data[k+2] / 255.0, 2.2);
			float r = in.data[k+0] / 255.0;
			float g = in.data[k+1] / 255.0;
			float b = in.data[k+2] / 255.0;

			float lum = r*0.35 + g*0.5 + b*0.15;
			//float lum = b;
			//float dither = noise.data[(i*noise.channels) % noise.length] / 255.0;
			//float dither = pow(noise.data[m] / 255.0, 2.2);
			float dither = ((noise.data[m] / 255.0) * 1.0 - 0.5);
			//float dither = 0;
			//float dither = rand()/(float)RAND_MAX - 0.5;

			//float srgb = pow(lum + dither, 1.0/2.2);
			//float srgb = pow(lum + dither, 2.2);
			float srgb = lum + dither;
			unsigned pxval = srgb < 0.5;

			buffer |= pxval << (7 - bits);
			//buffer |= pxval << bits;
			bits++;

			//fprintf(out, "%d ", srgb < 0.5);
			if (bits == 8) {
				fputc(buffer, out);
				bits = 0;
				buffer = 0;
			}
		}
	}

	if (bits != 0) {
		fputc(buffer, out);
	}

	fclose(out);

	return 0;
}
