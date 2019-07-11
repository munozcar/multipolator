#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "multipolator.h"

int main() {

  int file_descriptor;
  size_t length;
  double *grid;                                 // Grid will have values of type double.
  struct stat sb;

  // Open file with data.
  char *file_dir = "ATMO_GRID.txt";
  file_descriptor = open(file_dir, O_RDONLY);

  // Calculate size of grid & check for reading error.
  if (fstat(file_descriptor, &sb) == -1) {exit(EXIT_FAILURE);}
  length = sb.st_size; // Length, in bytes.

  // Map file into RAM
  grid = mmap(NULL, length, PROT_READ, MAP_PRIVATE, file_descriptor, 0);
  if (grid == MAP_FAILED) {exit(EXIT_FAILURE);}


  int iteration = 0;
  double interpolation_params[6] = {600, 6, 0.6, 0.6, 5, 0.6};
  double interpolated_model[5000];

  int param_N = grid[0];      // Number of parameters (dimensions) in the grid
  int points_N = grid[1];     // (Max) Number of model points per parameter


  while (iteration<1000){
    multipolator(grid, interpolation_params, interpolated_model, param_N, points_N);
    for (int i=0; i<param_N; i++){
      interpolation_params[i]*=1.00001;
    }
    for(int i=0; i<points_N; i++){
      fprintf(stdout, "DONE.");
    }
    iteration++;
  }

  // Remove RAM mapping, close file. Fin.
  munmap(grid, length);
  close(file_descriptor);
  return(EXIT_SUCCESS);
}
