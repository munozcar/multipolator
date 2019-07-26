#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include "multipolator.h"

int main() {

  int file_descriptor;
  size_t length;
  double *grid;                                 // Grid will have values of type double.
  struct stat sb;

  // Open file with data.
  char *file_dir = "ATMO_GRID_NOMALIZED.txt";
  file_descriptor = open(file_dir, O_RDONLY);

  // Calculate size of grid & check for reading error.
  if (fstat(file_descriptor, &sb) == -1) {printf("error"); exit(EXIT_FAILURE);}
  length = sb.st_size; // Length, in bytes.

  // Map file into RAM
  grid = mmap(NULL, length, PROT_READ, MAP_PRIVATE, file_descriptor, 0);
  if (grid == MAP_FAILED) {printf("error"); exit(EXIT_FAILURE);}

  double interpolation_params[6];

  interpolation_params[0] = 0.6818181818181818;
  interpolation_params[1] = 0.5555555555555556;
  interpolation_params[2] = 0.6391304347826088;
  interpolation_params[3] = 0.36153846153846153;
  interpolation_params[4] = 0.13557779799818018;
  interpolation_params[5] = 0.2;

  double interpolated_model[5000];

  //clock_t t;
  //t = clock();
  //for (int i=0; i<1000; i++){
  multipolator(grid, interpolation_params, interpolated_model);
  //}

  //t = clock() - t;
  //double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

  //printf("50000 mutlipolator() calls took %f seconds to execute \n", time_taken);
  printf("done");
  FILE *f;
  f = fopen("./degoutput.txt", "w");
  for (int i=0; i<5000; i++){
    fprintf(f, "%.12lf\n", interpolated_model[i]);
  }
  // Remove RAM mapping, close file. Fin.
  fclose(f);
  munmap(grid, length);
  close(file_descriptor);
  return(EXIT_SUCCESS);
}
