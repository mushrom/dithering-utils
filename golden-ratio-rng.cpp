#include <cstdio>
#include <cmath>
#include <ctime>

#define NUM_ITERS 20

int choice(float pick) {
	if (pick < 0.4) return 0;
	if (pick < 0.6) return 1;
	if (pick < 0.9) return 2;
	else            return 3;
}

int main(void) {
	srand(time(NULL));
	float seed = rand() / (float)RAND_MAX;
	float grat = (1 + sqrt(5)) / 2;

	unsigned counts[4]; 

	for (size_t i = 0; i < 4; i++) {
		counts[i] = 0;
	}

	for (size_t i = 0; i < NUM_ITERS; i++) {
		int c = choice(seed);

		printf("%3lu: %.4f: %8d\n", i, seed, c);
		counts[c]++;
		seed += grat;
		seed -= floor(seed);
	}

	for (size_t i = 0; i < 4; i++) {
		printf("%3lu: %.2f\n", i, counts[i] / (float)NUM_ITERS);
	}

	return 0;
}
