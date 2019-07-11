// ------------------------------------------------------------------------------------------------
// MULTIPOLATOR: Multidimensional Interpolator
// Written and updated by Carlos E. Munoz-Romero, 2019
// ------------------------------------------------------------------------------------------------

#define MAX_PARAMS 10
#define MAX_PARAMVALS 100

// FUNCTION PROTOTYPES ----------------------------------------------------------------------------

int factorial(int num);

int choose(int n, int k);

void swap(int *elements, int i, int j);

void reverse(int *elements, int i, int j);

void find_permutations(int *elements, int last_index, int permutation_count, int permutation_indices[][last_index+1], int *pnum);

void multipolator(double *grid, double *interpolation_parameters, double *model);
