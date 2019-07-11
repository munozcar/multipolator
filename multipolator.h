// ------------------------------------------------------------------------------------------------
// MULTIPOLATOR: Multidimensional Interpolator
// Written and updated by Carlos E. Munoz-Romero, 2019
// ------------------------------------------------------------------------------------------------

#define MAX_PARAMS 10
#define MAX_PARAMVALS 100
#define IGNORE_PRINTF 1

#ifdef IGNORE_PRINTF
#define printf(fmt, ...) (0)
#endif
// FUNCTION PROTOTYPES ----------------------------------------------------------------------------

int factorial(int num);

int choose(int n, int k);

void swap(int *elements, int i, int j);

void reverse(int *elements, int i, int j);

void find_permutations(int *elements, int last_index, int permutation_count, int permutation_indices[][last_index+1], int *pnum);

void multipolator(double *grid, double *parameters, double *model, int param_N, int points_N);
