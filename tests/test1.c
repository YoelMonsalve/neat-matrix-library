#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "lib/nml.h"
#include "lib/nml_def.h"
#include "nml_test.h"

#define TEST_TOLERANCE 0.0001

int main(int argc, char *argv[]) {
    
    // zeroes matrix
    nml_mat *A = nml_mat_new(3,4);
    nml_mat_printf(A, "%10.4lf  ");

    // eye matrix
    nml_mat *B = nml_mat_eye(4);
    printf("B: (%d) x (%d)\n", 
        nml_mat_dim(B, 1),
        nml_mat_dim(B, 2)
        );
    printf("B: (%d) x (%d)\n", 
        nml_mat_num_rows(B),
        nml_mat_num_cols(B)
        );
    nml_mat_printf(B, "%10.4lf  ");

    // random matrix
    nml_mat *C = nml_mat_rnd(3, 4, -1.0, +1.0);
    nml_mat_printf(C, "%10.4lf  ");    

    nml_mat_free(A);
    nml_mat_free(B);
    nml_mat_free(C);

    /*FILE *t_file = fopen(argv[1], "r");
    if (NULL == t_file) {
        fprintf(stderr, "Cannot open input file. %s", argv[1]);
        return -1;
    }
    int k, num_cases;
    fscanf(t_file, "%d", &num_cases);
    NML_TEST_HEADER(num_cases, argv[1], argv[0]);
    for(k = 0; k < num_cases; k++) {  
        nml_mat *A = nml_mat_fromfilef(t_file);
        nml_mat_lup *A_LUP = nml_mat_lup_solve(A);
        nml_mat *B = nml_mat_fromfilef(t_file);
        nml_mat *x_expected = nml_mat_fromfilef(t_file);
        nml_mat *x_computed = nml_ls_solve(A_LUP, B);

        if (nml_mat_eq(x_expected, x_computed, TEST_TOLERANCE)) {
            printf(GREEN "%d " RESET, k); 
        } else {
            printf(RED "%d" RESET, k);
            printf("A=\n");
            nml_mat_print(A);
            printf("A_LUP=\n");
            nml_mat_lup_print(A_LUP);
            printf("B=\n");
            nml_mat_print(B);
            printf("x_expected=\n");
            nml_mat_print(x_expected);
            printf("x_computed=\n");
            nml_mat_print(x_computed);
        }

        nml_mat_free(A);
        nml_mat_lup_free(A_LUP);
        nml_mat_free(B);
        nml_mat_free(x_computed);
        nml_mat_free(x_expected);
    }
    printf("\n");
    fclose(t_file);*/
    
    return 0;
}