#include <cstdio>
#include <cmath>
#include <ctime>

#include <vector>
#include <array>

constexpr size_t SIZE  = 32;
constexpr size_t TILES = 1;
constexpr float  SIGMA = 1.2f;
constexpr int    BOX   = std::ceil(2*SIGMA*SIGMA + 1);
constexpr float  SAMPLE_DENSITY = 0.25;

typedef std::array<std::array<float, 2*BOX+1>, 2*BOX+1> kernel;

struct vec2 {
	union {
		struct { float x, y; };
		float data[2];
	};
} __attribute__((packed));

struct ivec2 {
	union {
		struct { int x, y; };
		int data[2];
	};
} __attribute__((packed));

template <int B = BOX>
struct make_kernel {
	constexpr make_kernel() : data() {
		for (int y = -BOX; y < BOX; y++) {
			for (int x = -BOX; x < BOX; x++) {
				float r = sqrt(x*x + y*y);
				float fact = pow(M_E, -((r*r)/(2*(SIGMA*SIGMA))));

				data[y+BOX][x+BOX] = fact;
			}
		}
	}

	kernel data;
};

constexpr kernel guassian_kernel = make_kernel<>().data;

float random_float() {
	return rand() / (float)RAND_MAX;
}

vec2 random_vec2() {
	return (vec2) {
		.x = random_float(),
		.y = random_float(),
	};
}

float random_grd(float& seed) {
	constexpr float grat = (1 + std::sqrt(5)) / 2;

	seed += grat;
	seed -= floor(seed);

	return seed;
}

ivec2 random_offset() {
	//static float seed_x = random_float();
	//static float seed_y = random_float();
	static float seed_x = 0.001;
	static float seed_y = 0.5;

	//float v = random_grd();
	//float vx = random_grd(seed_x);
	//float vy = random_grd(seed_y);
	float vx = random_float();
	float vy = random_float();
	//fprintf(stderr, "%.2f, %.2f\n", vx, vy);

	return (ivec2) {
		/*
		.x = 0,
		.y = 0,
		*/
		.x = (int)(vx * SIZE),
		.y = (int)(vy * SIZE),
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

//#define SAMPLE_CANDIDATES 8
//#define NUM_SAMPLES (4*SIZE)
constexpr size_t SAMPLE_CANDIDATES = 8;
constexpr size_t NUM_SAMPLES = SAMPLE_DENSITY * (SIZE*SIZE);
vec2 samples[NUM_SAMPLES];
unsigned num_samples = 0;
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

std::vector<bool> gen_input_image_white_noise() {
	std::vector<bool> ret;
	unsigned ones = 0;

	for (size_t y = 0; y < SIZE; y++) {
		for (size_t x = 0; x < SIZE; x++) {
			bool val = random_float() < SAMPLE_DENSITY;
			ret.push_back(val);
			ones += val == true;
		}
	}

	// if there are more ones than zeros, flip so there are less
	if (ones > SIZE*SIZE/2) {
		fprintf(stderr, "flipping ones, %u %lu", ones, SIZE*SIZE);
		for (size_t i = 0; i < ret.size(); i++) {
			ret[i] = !ret[i];
		}
	}

	return ret;
}

std::vector<bool> gen_input_image_grd() {
	std::vector<bool> ret;

	ret.reserve(SIZE*SIZE);
	ret.resize(SIZE*SIZE);

	for (size_t i = 0; i < NUM_SAMPLES; i++) {
		ivec2 v = random_offset();
		ret[v.y*SIZE + v.x] = true;
	}

	return ret;
}


std::vector<bool> gen_input_image_mitchells() {
	std::vector<bool> ret;
	unsigned ones = 0;

	num_samples = 0;
	/*
	for (size_t y = 0; y < SIZE; y++) {
		for (size_t x = 0; x < SIZE; x++) {
			bool val = random_float() < SAMPLE_DENSITY;
			ret.push_back(val);
			ones += val == true;
		}
	}
	*/

	ret.reserve(SIZE*SIZE);
	ret.resize(SIZE*SIZE);

	for (size_t i = 0; i < NUM_SAMPLES; i++) {
		vec2 v = next_sample_vec2();

		unsigned x = v.x * SIZE;
		unsigned y = v.y * SIZE;

		ret[y*SIZE + x] = true;
		ones++;
	}

	/*
	// if there are more ones than zeros, flip so there are less
	if (ones > SIZE*SIZE/2) {
		fprintf(stderr, "flipping ones, %u %lu", ones, SIZE*SIZE);
		for (size_t i = 0; i < ret.size(); i++) {
			ret[i] = !ret[i];
		}
	}
	*/

	return ret;
}

int mod_size(int n, int inc) {
	if (n + inc < 0) {
		return SIZE + (n + inc);
	} else {
		return (n + inc) % SIZE;
	}
}

ivec2 mod_ivec2(ivec2 pos, ivec2 inc) {
	return (ivec2) {
		.x = mod_size(pos.x, inc.x),
		.y = mod_size(pos.y, inc.y),
	};
}

ivec2 mod_ivec2(ivec2 pos) {
	return (ivec2) {
		.x = (int)pos.x % (int)SIZE,
		.y = (int)pos.y % (int)SIZE,
	};
}

size_t index_ivec2(ivec2 idx) {
	return idx.y*SIZE + idx.x;
}

float convolve(std::vector<bool>& input, ivec2 pos) {
	float ret = 0;

	for (int k = -BOX; k < BOX; k++) {
		for (int m = -BOX; m < BOX; m++) {
			ivec2 idx = mod_ivec2(pos, {k, m});
			/*
			float r = sqrt(k*k + m*m);
			float fact = pow(M_E, -((r*r)/(2*(SIGMA*SIGMA))));
			*/
			float fact = guassian_kernel[m+BOX][k+BOX];

			ret += input[index_ivec2(idx)] * fact;
		}
	}

	return ret;
}

ivec2 find_largest_void(std::vector<bool>& input) {
	ivec2 ret;
	ivec2 off = random_offset();
	float cur = 99999;

	for (int y = 0; y < (int)SIZE; y++) {
		for (int x = 0; x < (int)SIZE; x++) {
			ivec2 offidx = mod_ivec2({.x = off.x + x, .y = off.y + y});
			bool v = input[index_ivec2(offidx)];

			if (v == 0) {
				float phi = convolve(input, offidx);

				if (phi < cur) {
					ret = offidx;
					cur = phi;
				}
			}
		}
	}

	return ret;
}

ivec2 find_largest_cluster(std::vector<bool>& input) {
	ivec2 ret;
	ivec2 off = random_offset();
	float cur = -1;

	for (int y = 0; y < (int)SIZE; y++) {
		for (int x = 0; x < (int)SIZE; x++) {
			ivec2 offidx = mod_ivec2({.x = off.x + x, .y = off.y + y});
			bool v = input[index_ivec2(offidx)];

			if (v == 1) {
				float phi = convolve(input, offidx);

				if (phi > cur) {
					ret = offidx;
					cur = phi;
				}
			}
		}
	}

	return ret;
}

bool operator==(const ivec2& a, const ivec2& b) {
	return a.x == b.x && a.y == b.y;
}

std::vector<bool> gen_initial_binary(std::vector<bool>& input) {
	std::vector<bool> ret;
	ivec2 last_removed = {-1, -1};

	ret.insert(ret.end(), input.begin(), input.end());

	for (;;) {
		ivec2 largest_void    = find_largest_void(ret);
		ivec2 largest_cluster = find_largest_cluster(ret);

		if (largest_void == last_removed) {
			ret[index_ivec2(last_removed)] = 1;
			//fprintf(stderr, "last removed");
			return ret;

		} else {
			/*
			fprintf(stderr, "swapping: (%d, %d) and (%d, %d)\n",
			        largest_void.x, largest_void.y, largest_cluster.x, largest_cluster.y);
					*/

			ret[index_ivec2(largest_void)]    = 1;
			ret[index_ivec2(largest_cluster)] = 0;
			last_removed = largest_cluster;
		}
	}

	return ret;
}

std::vector<int> gen_dither_matrix(std::vector<bool>& input) {
	std::vector<bool> prototype;
	std::vector<int> dither;

	dither.resize(input.size());

	unsigned ones = 0;
	int rank = 0;

	for (size_t i = 0; i < input.size(); i++) {
		ones += input[i];
	}

	// phase 1
	prototype = input;
	rank = ones - 1;

	while (rank >= 0) {
		ivec2 largest_cluster = find_largest_cluster(prototype);
		size_t idx = index_ivec2(largest_cluster);

		prototype[idx] = 0;
		dither[idx] = rank;
		rank--;
	}

	// phase 2
	// TODO: wait, why do we need phase 3? Maybe I'm misunderstanding this (likely),
	//       but it should behave exactly the same as finding the largest void no?
	prototype = input;
	rank = ones;

	while (rank < (int)(SIZE*SIZE)) {
		ivec2 largest_void = find_largest_void(prototype);
		size_t idx = index_ivec2(largest_void);

		prototype[idx] = 1;
		dither[idx] = rank;
		rank++;
	}

	return dither;
}

void dump_bool(std::vector<bool>& image) {
	printf("P3\n%lu %lu\n255\n", TILES*SIZE, TILES*SIZE);

	for (size_t y = 0; y < TILES*SIZE; y++) {
		for (size_t x = 0; x < TILES*SIZE; x++) {
			if (image[(y%SIZE)*SIZE + (x % SIZE)]) {
				puts("255 255 255");
			} else {
				puts("0 0 0");
			}
		}
	}
}

unsigned quantize(int rank) {
	return std::min(255.f, (256.f / (SIZE*SIZE)) * (rank + 0.5f));
	//return rank % 256;
}

void dump_dithered_scalar(std::vector<int>& image) {
	printf("P6\n%lu %lu\n255\n", TILES*SIZE, TILES*SIZE);

	for (size_t y = 0; y < TILES*SIZE; y++) {
		for (size_t x = 0; x < TILES*SIZE; x++) {
			int v = image[(y%SIZE)*SIZE + (x % SIZE)];
			unsigned px = quantize(v);

			//printf("%u %u %u\n", px, px, px);
			putchar(px);
			putchar(px);
			putchar(px);
		}
	}
}

void dump_dithered_rgb(std::vector<int>& a,
                       std::vector<int>& b,
                       std::vector<int>& c)
{
	printf("P3\n%lu %lu\n255\n", TILES*SIZE, TILES*SIZE);

	for (size_t y = 0; y < TILES*SIZE; y++) {
		for (size_t x = 0; x < TILES*SIZE; x++) {
			int rv = a[(y%SIZE)*SIZE + (x % SIZE)];
			int gv = b[(y%SIZE)*SIZE + (x % SIZE)];
			int bv = c[(y%SIZE)*SIZE + (x % SIZE)];
			unsigned r = quantize(rv);
			unsigned g = quantize(gv);
			unsigned b = quantize(bv);

			printf("%u %u %u\n", r, g, b);
		}
	}
}

void dump_debug_info(std::vector<int>& image) {
	size_t counts[256];

	for (size_t i = 0; i < 256; i++) {
		counts[i] = 0;
	}

	for (size_t y = 0; y < SIZE; y++) {
		for (size_t x = 0; x < SIZE; x++) {
			int v = image[y*SIZE + x];
			unsigned px = quantize(v);

			counts[px]++;
		}
	}

	for (size_t i = 0; i < 256; i++) {
		fprintf(stderr, "%02lx: %lu (%.2f%%)   ", i, counts[i], 100.f*counts[i]/(float)(SIZE*SIZE));

		if ((i + 1) % 6 == 0) {
			fprintf(stderr, "\n");
		}
	}

	fprintf(stderr, "\n");
	fprintf(stderr, "box size: %d\n", BOX);
}

std::vector<int> gen_dither_channel() {
	//auto input   = gen_input_image_white_noise();
	//auto input   = gen_input_image_grd();
	auto input   = gen_input_image_mitchells();
	//dump_bool(input);
	//return {};
	auto initial = gen_initial_binary(input);
	auto dither  = gen_dither_matrix(initial);

	return dither;
}

int main(void) {
	srand(time(NULL));

	/*
	auto input   = gen_input_image_white_noise();
	auto initial = gen_initial_binary(input);
	auto dither  = gen_dither_matrix(initial);
	*/

	/*
	auto dither_r = gen_dither_channel();
	auto dither_g = gen_dither_channel();
	auto dither_b = gen_dither_channel();

	dump_dithered_rgb(dither_r, dither_g, dither_b);
	*/

	auto dither = gen_dither_channel();
	//return 0;

	//dump_bool(input);
	//dump_bool(initial);
	dump_dithered_scalar(dither);
	dump_debug_info(dither);

	return 0;
}
