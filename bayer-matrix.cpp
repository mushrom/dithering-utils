#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdint>

#include <vector>
#include <array>

constexpr size_t SIZE = 8;
constexpr size_t TILES = 1;
constexpr size_t MULT = 256 / (SIZE * SIZE);

void dump_dithered_scalar(std::vector<int>& image) {
	printf("P6\n%d %d\n255\n", TILES*SIZE, TILES*SIZE);

	for (size_t y = 0; y < TILES*SIZE; y++) {
		for (size_t x = 0; x < TILES*SIZE; x++) {
			int px = image[(y%SIZE)*SIZE + (x % SIZE)];

			//printf("%u %u %u\n", px, px, px);
			putchar(px);
			putchar(px);
			putchar(px);
		}
	}
}

unsigned reverse(unsigned num) {
	unsigned ret = 0;

	for (size_t i = 0; i < 8; i++) {
		unsigned b = !!(num & (1 << i));
		ret |= b << (7 - i);
	}

	return ret;
}

void fill_matrix(float matrix[SIZE][SIZE], unsigned n, unsigned px, unsigned py) {
	if (n == 1)
		return;

	fill_matrix(matrix, n/2, px,       py);
	fill_matrix(matrix, n/2, px + n/2, py);
	fill_matrix(matrix, n/2, px,       py + n/2);
	fill_matrix(matrix, n/2, px + n/2, py + n/2);

	for (size_t y = py; y < py+n; y++) {
		for (size_t x = px; x < px+n; x++) {
			matrix[y][x] *= n*n;

			     if (x <  px+n/2 && y <  py+n/2) matrix[y][x] += 0.f;
			else if (x >= px+n/2 && y <  py+n/2) matrix[y][x] += 2.f;
			else if (x <  px+n/2 && y >= py+n/2) matrix[y][x] += 3.f;
			else if (x >= px+n/2 && y >= py+n/2) matrix[y][x] += 1.f;
			else fprintf(stderr, "wut?");

			matrix[y][x] *= 1.f/(n*n);
		}
	}
}

void dump_matrix(float matrix[SIZE][SIZE]) {
	for (size_t y = 0; y < SIZE; y++) {
		for (size_t x = 0; x < SIZE; x++) {
			fprintf(stderr, "%3.f ", (SIZE*SIZE)*matrix[y][x]);
		}

		fprintf(stderr, "\n");
	}
}

int main(void) {
	float matrix[SIZE][SIZE];
	for (size_t y = 0; y < SIZE; y++) {
		for (size_t x = 0; x < SIZE; x++) {
			matrix[y][x] = 0;
		}
	}

	fill_matrix(matrix, SIZE, 0, 0);
	dump_matrix(matrix);

	printf("P6\n%d %d\n255\n", SIZE, SIZE);

	for (size_t y = 0; y < SIZE; y++) {
		for (size_t x = 0; x < SIZE; x++) {
			int px = matrix[y][x] * 255;

			putchar(px);
			putchar(px);
			putchar(px);
		}
	}

	return 0;
}
