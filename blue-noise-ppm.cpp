#include <cstdio>
#include <cmath>
#include <ctime>

struct vec2 {
	union {
		struct { float x, y; };
		float data[2];
	};
} __attribute__((packed));

#define SIZE 64
#define TILES 8
#define SAMPLE_CANDIDATES 8
//#define NUM_SAMPLES (SIZE*SIZE)
#define NUM_SAMPLES (4*SIZE)

unsigned num_samples = 0;
vec2 samples[NUM_SAMPLES];

float random_float() {
	return rand() / (float)RAND_MAX;
}

vec2 random_vec2() {
	return (vec2) {
		.x = random_float(),
		.y = random_float(),
	};
}

float vec2_dist(vec2& a, vec2& b) {
	float x = fabs(a.x - b.x);
	float y = fabs(a.y - b.y);

	if (x > 0.5) x = 1.0 - x;
	if (y > 0.5) y = 1.0 - y;

	return sqrt(x*x + y*y);
	//return fabs(x) + fabs(y) + fabs(z);
}

vec2 next_sample_vec2() {
	vec2 cur = random_vec2();
	float dist = 0;

	for (size_t i = 0; i < SAMPLE_CANDIDATES*num_samples + 1; i++) {
		float avg_dist = 99999;
		vec2 temp = random_vec2();

		for (size_t k = 0; k < num_samples; k++) {
			//avg_dist += vec2_dist(temp, samples[k]);
			//avg_dist = std::min(avg_dist, vec2_dist(temp, samples[k]));
			avg_dist = std::min(avg_dist, vec2_dist(temp, samples[k]));
			//avg_dist = std::min(avg_dist, fabsf(temp.x - samples[k].x));
		}

		if (avg_dist > dist) {
			cur  = temp;
			dist = avg_dist;
		}
	}

	samples[num_samples] = cur;
	num_samples++;

	return cur;
}

bool image[SIZE * SIZE];

int main(void) {
	srand(time(NULL));

	printf("P3\n%d %d\n255\n", TILES*SIZE, TILES*SIZE);

	for (size_t i = 0; i < NUM_SAMPLES; i++) {
		vec2 v = next_sample_vec2();

		unsigned x = v.x * SIZE;
		unsigned y = v.y * SIZE;

		image[y*SIZE + x] = true;
	}

	for (size_t y = 0; y < TILES*SIZE; y++) {
		for (size_t x = 0; x < TILES*SIZE; x++) {
			/*
			vec2 n = next_sample_vec2();
			n.x = floor(n.x * 255);
			//n.y = floor(n.y * 255);

			printf("%d %d %d\n", (unsigned)n.x, (unsigned)n.x, (unsigned)n.x);
			*/

			if (image[(y%SIZE)*SIZE + (x % SIZE)]) {
				puts("255 255 255");
			} else {
				puts("0 0 0");
			}

			//printf("%d %d %d\n", (int)n, (int)n, (int)n);
		}
	}

	return 0;
}
