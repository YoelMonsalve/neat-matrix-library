#include <stdlib.h>
#include <stdio.h>

#include "nml.h"

int main(int argc, char *argv[]) {

  double v1[9] = {
    0.0, 1.0, 2.0,
    1.0, 2.0, 1.0,
    2.0, 7.0, 8.0
  };

  nml_mat *m1 = nml_mat_from(3, 3, 9, v1);
  printf("\nm1=\n");
  nml_mat_print(m1);
  nml_mat *refm1 = nml_mat_ref(m1);
  printf("\nrefm1=\n");
  nml_mat_print(refm1);

  nml_mat_free(m1);
  nml_mat_free(refm1);
  return 0;
}
