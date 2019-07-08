/* ----------------------------------------------------------------------
multipolator: An n-dimensional grid interpolator.
Written and updated by Carlos E. Munoz-Romero, Geronimo L. Villanueva.[1]
[1] Goddard Center for Astrobiology, NASA GSFC, Greenbelt, MD.
Last edited: July 2019

multipolator finds nearest grid points to a given input and approximates the function through a
weighted average of all possible permutations. Weights are normalized using a partition function.
*/

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define MAX_PARAMVALS 100
#define MAX_PARAMS 10

// AUXILIARY PROCEDURES
//-------------------------------------------------------------------------------------------------
int factorial(int num){
  int fac = num;
  if (num==0){
    return 1;
  }
  for (int i=num; i>1; i--){
    fac*=(i-1);
  }
  return fac;
}

int choose(int n, int k){
  // n choose k algorithm
  int c = factorial(n) / (factorial(k) * factorial(n-k));
  return c;
}

void swap(int *elements, int i, int j){
  int temp = 0;
  temp = elements[i];
  elements[i] = elements[j];
  elements[j] = temp;
}

void reverse(int *elements, int i, int j){
  int end = floor((j-i+1)/2);
  for (int offset=0; offset<end; offset++){
    swap(elements, i+offset, j-offset);
  }
}

void find_permutations(int *elements, int last_index, int permutation_count, int permutation_indices[][last_index+1], int *pnum){

  if (last_index < 1){
    return;
  }
  int i = last_index - 1;
  int j = 0;
  while (i>=0 && !(elements[i]<elements[i+1])){
    i--;
  }
  if (i<0){
    reverse(elements, 0, last_index);
  }
  else {
    j = last_index;
    while(j>(i+1) && !(elements[j] > elements[i])){
      j--;
    }
    swap(elements, i, j);
    reverse(elements, i+1, last_index);
  }
  printf("CASE %d: ",pnum[0]);
  for (int k=0; k<last_index+1; k++){
    permutation_indices[pnum[0]][k] = elements[k];
    printf("%d ", permutation_indices[pnum[0]][k]);
  }
  pnum[0] = pnum[0] + 1;
  printf("\n");
  if (permutation_count>0){
    permutation_count--;
    find_permutations(elements, last_index, permutation_count, permutation_indices, pnum);
  }

}
// -----------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  if (argc != 5) {printf("\nWrong number of arguments! Exit.\n\n"); return EXIT_FAILURE;}

  int i = 0;
  int file_descriptor;
  size_t length;
  double *grid;                                 // Grid will have values of type double.
  int param_N = strtol(argv[1], NULL, 10);      // Number of parameters (dimensions) in the grid
  int paramvals_N = strtol(argv[2], NULL, 10);  // (Max) Number of parameter values in the grid
  int points_N = strtol(argv[3], NULL, 10);     // (Max) Number of model points per parameter
  struct stat sb;

  if (paramvals_N > MAX_PARAMVALS || param_N > MAX_PARAMS)
    {printf("\nParameter space too large! Exit.\n\n"); return EXIT_FAILURE;}

  double parameters[MAX_PARAMS][MAX_PARAMVALS] = {0};
  double interpolation_params[param_N];
  int interpolation_indices[MAX_PARAMS][2] = {0}; // Array for indices of 'left' and 'right' closest points.

  int permutations = (int)pow(2,param_N);  // Total number of permutations 2^(number of parameters)
  int multiplicatives[param_N];
  int param_shape[param_N];                // Exactly how many values for each parameter we have
  double models[permutations][points_N];   // Array to save models for interpolation

  // Open file with data.
  file_descriptor = open(argv[4], O_RDONLY);

  // Calculate size of grid & check for reading error.
  if (fstat(file_descriptor, &sb) == -1) {exit(EXIT_FAILURE);}
  length = sb.st_size; // Length, in bytes.

  // Map file into RAM
  grid = mmap(NULL, length, PROT_READ, MAP_PRIVATE, file_descriptor, 0);
  if (grid == MAP_FAILED) {exit(EXIT_FAILURE);}

  // Read the parameter space, shifting to the next parameter when a nan is encountered.
  // Skip first row, as it contains the covariates.

  for (int i=0; i < param_N; i++) {
    param_shape[i] = 0;
    for (int j=0; j < paramvals_N; j++){
      parameters[i][j] = grid[points_N*(i+1)+j];
      if (parameters[i][j] == parameters[i][j]) { // If not nan, we have not read all values of this parameter
        param_shape[i]++;
      }
    }
  }

  // Print parameter space
  printf("\n-----------------------------------------------------------------\n");
  printf("PARAMETER SPACE\n\nFOUND:\n\n");
  for (int j=0; j < param_N; j++){
     printf("%2d VALUES\t", param_shape[j]);
   }
  printf("\n");
  for (int i = 0; i < paramvals_N; ++i) {
    for (int j=0; j < param_N; j++){
      if (parameters[j][i] == parameters[j][i]) {
        printf("%8.2lf \t", parameters[j][i]);
      }
      else {
        printf("\t\t");
      }
     }
     printf("\n");
   }

   // We will interpolate to these parameters
   printf("\n-----------------------------------------------------------------\n");
   printf( "INTERPOLATION PARAMETERS:\n");

   for (int i=0; i < param_N; i++) {
     scanf("%lf", &interpolation_params[i]);
   }
   // Confirm parameters to interpolate
   printf("\n");
   printf( "WILL INTERPOLATE TO:\n");
   for (int i=0; i < param_N; i++) {
     printf("%8.2lf\t", interpolation_params[i]);
   }
   printf("\n-----------------------------------------------------------------\n");

  // Now we perform the interpolation

  // Find INDICES of the two closest points to each interpolation parameter
  printf( "PERMUTATION MODULE\n\nFINDING NEAREST TWO GRID POINTS TO EACH PARAMETER:\n");
  for (int i=0; i<param_N; i++) {
    while (parameters[i][interpolation_indices[i][0]] <= interpolation_params[i]) {
      interpolation_indices[i][0]++;
    }
    if (interpolation_indices[i][0]>0) {
      interpolation_indices[i][0]=interpolation_indices[i][0]-1;
    }
    interpolation_indices[i][1] = interpolation_indices[i][0]+1;
    // Check if input is out of grid range
    if (parameters[i][interpolation_indices[i][1]] != parameters[i][interpolation_indices[i][1]] ||
        interpolation_params[i] < parameters[i][0])
        {
          printf("\nParameter out of range! %8.2lf \n Exit.\n\n", interpolation_params[i]);
          return EXIT_FAILURE;
        }
  }
  for (int i=0; i < param_N; i++) {
    printf("%8.2lf\t%8.2lf\n", parameters[i][interpolation_indices[i][0]],
                               parameters[i][interpolation_indices[i][1]]);
  }

  // Next, we find and read all models corresponding to the permutations
  printf( "\nFINDING INDEX MULTIPLICATIVES:\n");

  // This step allows us to find the indices of the models of interest
  for (int i=0; i<param_N; i++){
    multiplicatives[i] = 1;
    for (int j=0; j<i; j++){
      multiplicatives[i]*=param_shape[j];
    }
    printf("%d\n", multiplicatives[i]);
  }

  // Calculate all posible permutations
  printf( "\nTHERE WILL BE (%d) MODEL PERMUTATIONS IN TOTAL.\n\nPASCAL MATRIX (0 AND 1 FOR LEFT AND RIGHT ENDPOINTS):\n\n", permutations);

  int permutation_indices[permutations][param_N];   // Will store all possible combinations here
  int permutation_count[param_N+1];
  int layer[param_N];                               // A single layer
  int layer_matrix[param_N+1][param_N];               // Each row will have its own permutations.
  int temp, skip_index;
  int curr = 0;
  int pnum[1] = {0};


  for (int i=0; i<param_N; i++) {layer[i] = 0;}     // Populate base case with zeros.

  for (int i=0; i<param_N+1; i++) {                   // Populate matrix
    for (int j=0; j<param_N; j++) {
      if (i<=j) {layer_matrix[i][j] = 1;}
      else {layer_matrix[i][j] = 0;}
      printf("%d ", layer_matrix[i][j]);
    }
    permutation_count[i] = choose(param_N,param_N-i);
    printf(" PERMUTATIONS: %d\n", permutation_count[i]);
  }
  printf("\n");

  printf("GENERATING PERMUTATIONS:\n\n");
  for (int i=0; i<param_N+1; i++){
    for (int j=0; j<param_N; j++){
      layer[j] = layer_matrix[i][j];
    }
    find_permutations(layer, param_N-1, permutation_count[i]-1, permutation_indices, pnum);
  }

  printf("\n FOR EACH CASE, SAVING MODELS...\n\n");

  for(int i=0; i<permutations; i++) {
    skip_index = 0;
    for(int j=0; j<param_N; j++) { // Find model index
      skip_index+=multiplicatives[j]*interpolation_indices[j][permutation_indices[i][j]];
    }
    skip_index*=points_N;
    for(int k=0; k<points_N; k++) { // Save model
      models[i][k] = grid[points_N*(param_N+1) + skip_index + k];
    }
  }

  printf("\n SUCCESS\n\n");

  printf("\n-----------------------------------------------------------------\n");
  printf("\n INTERPOLATION MODULE\n\n");

  

  printf( "\nDONE\n");

  // Remove RAM mapping, close file. Fin.
  munmap(grid, length);
  close(file_descriptor);

  return EXIT_SUCCESS;
}
