/**
Copyright 2021 Andrei N. Ciobanu

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "nml.h"

/** errors/exceptions */
#define DEFAULT_VALUE 0.0

#define CANNOT_ADD "Cannot add two matrices with different dimensions.\n"

#define CANNOT_SUBSTRACT "Cannot substract two matrices with different dimensions.\n"

#define CANNOT_MULTIPLY \
  "Cannot multiply two matrices where \
  the number of columns of the first one \
  is different than the number of rows of the second one.\n" \

#define CANNOT_REMOVE_COLUMN "Cannot remove matrix column %d. The value should be less than %d.\n" 

#define CANNOT_REMOVE_ROW "Cannot remove matrix row %d. The value should be less than %d.\n" 

#define INVALID_ARGUMENT_INT "Invalid argument: %d\n"

#define INVALID_ROWS \
  "Cannot create matrix with 0 number of rows. Aborting.\n" \

#define INVALID_COLS \
    "Cannot create matrix with 0 number of cols. Aborting.\n" \

#define CANNOT_TRACE \
    "Cannot calculate trace. Matrix needs to be square.\n" \

#define CANNOT_CROUT \
    "Cannot apply crout algorithm. Matrix needs to be square.\n" \

#define CANNOT_SWAP_ROWS \
     "Cannot swap rows (%d, %d) because the matrix number of rows is %d.\n" \

#define CANNOT_SWAP_COLUMNS \
      "Cannot swap columns (%d, %d) because the matrix number or columns is %d.\n" \

#define CANNOT_ROW_MULTIPLY \
      "Cannot multiply row (%d), maximum number of rows is %d.\n" \

#define CANNOT_COL_MULTIPLY "Cannot multiply col (%d), maximum number of columns is %d.\n" 
  
#define CANNOT_ADD_TO_ROW \
      "Cannot add %2.2f x (row=%d) to row=%d. Total number of rows is: %d.\n" \

#define CANNOT_LU_MATRIX_SQUARE \
      "Canot LU. Matrix (%d, %d) needs to be square.\n"

#define CANNOT_LU_MATRIX_DEGENERATE \
      "Cannot LU. Matrix is degenerate or almost degenerate.\n"

#define INCONSISTENT_SYSTEM \
      "Cannot LU. Matrix is degenerate or almost degenerate.\n"

#define CANNOT_SOLVE_LIN_SYS_INVALID_B \
      "Cannot solve system. b[%d][%d] should have size b[%d][%d].\n" \

#define CANNOT_SET_DIAG \
      "Cannot set diag with value(=%2.2f). Matrix is not square.\n" \

#define CANNOT_CONCATENATE_H \
      "Cannot concatenate. Matrices have a different number of rows. Expected %d, found: %d.\n" \

#define CANNOT_CONCATENATE_V \
      "Cannot concatenate. Matrices have a different number of cols. Expected %d, found: %d.\n" \

#define CANNOT_GET_COLUMN \
      "Cannot get column (%d). The matrix has %d number of columns.\n" \

#define CANNOT_GET_ROW \
      "Cannot get row (%d). The matrix has %d number of rows.\n" \

#define INCONSISTENT_ARRAY \
      "Cannot found element %d in the array (NULL). Expected a total of : %d elements.\n"  \

#define INCONSISTENT_VARGS \
      "Cannot find element %d in the varargs. Expecteda total of : %d varargs.\n" \

#define CANNOT_REF_MATRIX_DEGENERATE \
      "Cannot compute REF. Matrix is degenerate or near degenerate.\n" \

#define CANNOT_OPEN_FILE "Cannot open file '%s'. Please check the path is correct and you have reading rights.\n"

#define INVALID_MATRIX_FILE \
      "Invalid matrix file: %s. Cannot read data.\n" \

#define VECTOR_J_DEGENERATE \
      "Vector on colum %d is generate or near degenerate. Cannot proceed further.\n"

#define CANNOT_QR_NON_SQUARE \
      "We cannot QA non-square matrix[%d, %d].\n"

#define CANNOT_COLUMN_L2NORM \
      "Cannot get column (%d). The matrix has %d numbers of columns.\n"

#define CANNOT_VECT_DOT_DIMENSIONS \
      "The two vectors have different dimensions: %d and %d.\n"

#define MATRIX_NOT_SQUARE \
      "Matrix (%d, %d) is not square.\n"
       

// *****************************************************************************
//
// Constructing and destroying a matrix struct
//
// *****************************************************************************

/**
 * Dynamically allocates a new matrix struct.
 * The matrix will be initially populated with zeroes.
 *
 * @param   num_rows  number of rows of the new matrix
 * @param   num_cols  number of cols of the new matrix
 * @return            a pointer to the new matrix struct
 */
nml_mat *nml_mat_new(unsigned int num_rows, unsigned int num_cols) {
  if (num_rows == 0) {
    NML_ERROR(INVALID_ROWS);
    return NULL;
  }
  if (num_cols == 0) {
    NML_ERROR(INVALID_COLS);
    return NULL;
  }
  nml_mat *m = calloc(1, sizeof(*m));
  NP_CHECK(m);    // asserting the data was allocated
  
  m->num_rows = num_rows;
  m->num_cols = num_cols;
  
  m->is_square = (num_rows == num_cols) ? 1 : 0;
  m->data = calloc(m->num_rows, sizeof(*m->data));
  NP_CHECK(m->data);
  int i;
  for(i = 0; i < m->num_rows; ++i) {
    m->data[i] = calloc(m->num_cols, sizeof(**m->data));
    NP_CHECK(m->data[i]);
  }

  // new form, to store data contiguously
  m->__data = calloc(m->num_rows * m->num_cols, sizeof(*m->__data));
  NP_CHECK(m->__data);    // asserting the data was allocated
  return m;
}

/**
 * Free a matrix structure.
 *
 * @param  matrix  the matrix to be freed.
 */
void nml_mat_free(nml_mat *matrix) {
  int i;
  for(i = 0; i < matrix->num_rows; ++i) {
    free(matrix->data[i]);
    // Yoel.-
    matrix->data[i] = NULL;
  }
  free(matrix->data);
  free(matrix);
  // Yoel.-  It is a good practice to set the unused pointers to NULL.
  matrix->data = NULL;
  matrix = NULL;
}

/**
 * Dynamically allocates a new square matrix struct.
 * This is similar to `nml_matrix_new()`, but the matrix is square 
 * (number of rows is equal than the number of cols)
 */
nml_mat *nml_mat_sqr(unsigned int size) {
  return nml_mat_new(size, size);
}

/**
 * Create a new square matrix, with random elements.
 * The elements will be double in the range [min, max)
 *
 * @param   size  the size of the matrix
 * @param   min   lower bound for the random numbers
 * @param   max   upper bound for the random numbers
 * @return        a pointer to the created matrix struct
 */
nml_mat *nml_mat_rnd(unsigned int num_rows, unsigned int num_cols, double min, double max) {
  nml_mat *r = nml_mat_new(num_rows, num_cols);
  int i, j;
  
  #if 1
  for(i = 0; i < num_rows; i++) {
    for(j = 0; j < num_cols; j++) {
      r->data[i][j] = nml_rand_interval(min, max);
    }
  }
  #endif

  // new form, based in contiguous allocation
  int k;
  for (k = 0; k < num_rows * num_cols; k++) {
    __NML_DATA(r)[k] = nml_rand_interval(min, max);
  }
  return r;
}

/**
 * Create a new square matrix, with random elements.
 * The elements will be double in the range [min, max)
 *
 * @param   size  the size of the matrix
 * @param   min   lower bound for the random numbers
 * @param   max   upper bound for the random numbers
 * @return        a pointer to the created matrix struct
 */
nml_mat *nml_mat_sqr_rnd(unsigned int size, double min, double max) {
  return nml_mat_rnd(size, size, min, max);
}

/**
 * Dynamically allocates a new identity matrix.
 * The resulting matrix will be sqare, with 1's in the main diagonal
 * and zero otherwise.
 *
 * @param   size  the size of the square matrix
 * @return        identity matrix
 */
nml_mat *nml_mat_eye(unsigned int size) {
  nml_mat *r = nml_mat_new(size, size);
  int i;
  for(i = 0; i < r->num_rows; i++) {
    r->data[i][i] = 1.0;
  
    // new form, based on contiguous allocation
    const int num_cols = __NML_COLS(r);
    __NML_ELEM(r, i, i, num_cols) = 1.0;
  }
  return r;
}

/**
 * Initialise a matrix by reading values from a vector.
 *
 * @param   num_rows  number of rows of the new matrix
 * @param   num_cols  number of columns of the new matrix
 * @param   n_vals    lenght of the vector of values to read from
 * @param   vals      vector of values to read from
 * @return            pointer to the new created matrix struct.
 */
nml_mat *nml_mat_from(unsigned int num_rows, unsigned int num_cols, unsigned int n_vals, double *vals) {
  int i, j, v_idx;

  nml_mat *m = nml_mat_new(num_rows, num_cols);
  
  #if 1
  // Yoel.-  m->num_rows is the same than num_rows, but using num_rows will avoid the
  //         overload by accessing m.
  //for(i = 0; i < m->num_rows; i++) {
  for(i = 0; i < num_rows; i++) {
    //for(j = 0; j < m->num_cols; j++) {
    for(j = 0; j < num_cols; j++) {
      v_idx = i * num_cols + j;
      m->data[i][j] = (v_idx < n_vals) ? vals[v_idx] : 0.0;
    }
  }
  #endif

  /* Delegating on the efficient standard C library, the task to deal with
   * memory operations */
  memcpy(__NML_DATA(m), vals, n_vals);
  return m;
}

/**
 * Initialise a matrix by copying from another one.
 * Dynamically allocates a new Matrix.
 *
 * @param   m  the other matrix
 * @return     pointer to the new created matrix struct.
 */
nml_mat *nml_mat_cp(nml_mat *m) {
  int i,j;
  const unsigned int num_rows = m->num_rows,
    num_cols = m->num_cols;

  #if 1
  nml_mat *r  = nml_mat_new(num_rows, num_cols);
  for(i = 0; i < num_rows; i++) {
    for(j = 0; j < num_cols; j++) {
      r->data[i][j] = m->data[i][j];
    }
  }
  #endif

  /* Delegating on the efficient standard C library, the task to deal with
   * memory operations */
  memcpy(__NML_DATA(r), __NML_DATA(m), num_rows * num_cols);
  return r;
}

nml_mat *nml_mat_fromfile(const char *file) {
  FILE *m_file = fopen(file, "r");
  if (NULL == m_file) {
    NML_FERROR(CANNOT_OPEN_FILE, file);
    return NULL;
  }
  nml_mat *r = nml_mat_fromfilef(m_file);
  fclose(m_file);
  return r;
}

nml_mat *nml_mat_fromfilef(FILE *f) {
  int i, j;
  unsigned int num_rows = 0, num_cols = 0;
  fscanf(f, "%d", &num_rows);
  fscanf(f, "%d", &num_cols);
  nml_mat *r = nml_mat_new(num_rows, num_cols);
  for(i = 0; i < num_rows; i++) {
    for(j = 0; j < num_cols; j++) {
      //fscanf(f, "%lf\t", &r->data[i][j]);
      // new form, based in contiguous allocation
      // this is the same than
      //   &(r->__data[ i*num_cols + j ])
      // or
      //   r->__data + (i*num_cols + j)
      fscanf(f, "%lf\t", __NML_DATA(r) + __NML_1D_INDEX(i, j, num_cols));
    }
  }
  return r;
}

void nml_mat_from_array(nml_mat *A, const double *v) {
  int i, j;
  unsigned int num_rows = __NML_ROWS(A), num_cols = __NML_COLS(A);
  for (i = 0; i < num_rows; i++) {
    for (j = 0; j < num_cols; j++) {
      __NML_ELEM2(A, i, j) = v[i*num_cols + j];
    }
  }
}

/**
 * Return the number of rows of matrix.
 *
 * @param   m  pointer to matrix struct
 * @return     number of rows
 */
inline int nml_mat_num_rows(nml_mat *A) { return __NML_ROWS(A); }

/**
 * Return the number of cols of matrix.
 *
 * @param   m  pointer to matrix struct
 * @return     number of cols
 */
inline int nml_mat_num_cols(nml_mat *A) { return __NML_COLS(A); }

/**
 * Number of elements in the i-th dimension
 *
 * - dim 1: number of rows
 * - dim 2: number of columns
 */
inline int nml_mat_dim(nml_mat *A, int dim) {
  if (dim == 1) {
    return __NML_ROWS(A);
  }
  else if (dim == 2) {
    return __NML_COLS(A);
  }
  else {
    NML_FERROR(INVALID_ARGUMENT_INT, dim);
    return -1;
  }
}

// *****************************************************************************
//
// Matrix Equality
//
// *****************************************************************************

// Checks if two matrices have the same dimensions
// 
// As this is a very short code, it worth to be implemented
// as an inline function (C99).
inline int nml_mat_eqdim(nml_mat *m1, nml_mat *m2) {
  return (m1->num_cols == m2->num_cols) &&
          (m1->num_rows == m2->num_rows);
}

// Checks if two matrices have the same dimensions, and the elements
// are all equal to each other with a given tolerance;
// For exact equality use tolerance = 0.0
int nml_mat_eq(nml_mat *m1, nml_mat *m2, double tolerance) {
  if (!nml_mat_eqdim(m1, m2)) {
    return FALSE;    //  will not cause overhead, as MACROS will be substituted at preprocessing time
  }
  int i, j;
  #if 1
  for(i = 0; i < m1->num_rows; i++) {
    for(j = 0; j < m1->num_cols; j++) {
      if (fabs(m1->data[i][j] - m2->data[i][j]) > tolerance) {
        return FALSE;
      }
    }
  }
  #endif

  // new form, based on contiguous allocation
  const int num_rows = m1->num_rows, num_cols = m1->num_cols;
  int k;
  for (k = 0; k < num_rows * num_cols; k++) {
    if (fabs(__NML_DATA(m1)[k] - __NML_DATA(m2)[k]) > tolerance) {
      return FALSE;
    }
  }
  return TRUE;
}

// *****************************************************************************
//
// Matrix printing
//
// *****************************************************************************

// Prints the matrix on the stdout
void nml_mat_print(nml_mat *matrix) {
  nml_mat_printf(matrix, "%lf\t\t");
}

// Prints the matrix on the stdout (with a custom formatting for elements)
void nml_mat_printf(nml_mat *matrix, const char *d_fmt) {
  const int num_rows = __NML_ROWS(matrix), num_cols = __NML_COLS(matrix);
  int i, j;

  fprintf(stdout, "\n");
  for(i = 0; i < num_rows; ++i) {
    for(j = 0; j < num_cols; ++j) {
      //fprintf(stdout, d_fmt, matrix->data[i][j]);

      // new form, based on contiguous allocation
      fprintf(stdout, d_fmt, __NML_ELEM(matrix, i, j, num_cols));      
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "\n");
}

// *****************************************************************************
//
// Accessing and modifying matrix elements
//
// *****************************************************************************
/**
 * Get the [i,j]-th element of m (zero-based index)
 * C99 allows the qualifier "inline" to suggest the compiler to insert the 
 * code function in-place, instead of calling it as a normal function.
 *
 * @param   matrix  the matrix (a pointer to)
 * @param   i       row index
 * @param   j       col index
 * @return          the [i,j]-th element
 */
inline double nml_mat_get(nml_mat *matrix, unsigned int i, unsigned int j) {
  // new form, based in contiguous allocation
  //return matrix->data[i][j];
  return __NML_ELEM(matrix, i, j, matrix->num_cols);
}

/**
 * Set the [i,j]-th element of m (zero-based index)
 *
 * C99 allows the qualifier "inline" to suggest the compiler to insert the 
 * code function in-place, instead of calling it as a normal function.
 * 
 * @param   matrix  the matrix (a pointer to)
 * @param   i       row index
 * @param   j       col index
 * @param   val     the value to what m[i,j] will be set
 */
inline void nml_mat_set(nml_mat *matrix, unsigned int i, unsigned int j, double val) {
  //matrix->data[i][j] = val;
  // new form, based in contiguous allocation
  __NML_ELEM(matrix, i, j, matrix->num_cols) = val;
}

// produces a column-matrix [Nx1] from the col-th column of m
nml_mat *nml_mat_col_get(nml_mat *m, unsigned int col) {
  const unsigned int num_cols = m->num_cols,
    num_rows = m->num_rows;

  if (col >= num_cols) {
    NML_FERROR(CANNOT_GET_COLUMN, col, num_cols);
    return NULL;
  }
  nml_mat *r = nml_mat_new(num_rows, 1);
  int i;
  for(i = 0; i < num_rows; i++) {
    //r->data[i][0] = m->data[i][col];
    // new form, based in contiguous allocation
    __NML_ELEM(r, i, 0, 1) = __NML_ELEM(m, i, col, num_cols);
  }
  return r;
}

// produces a row-matrix [1xN] from the row-th row of m
nml_mat *nml_mat_row_get(nml_mat *m, unsigned int row) {
  const unsigned int num_rows = m->num_rows,
    num_cols = m->num_cols;

  if (row >= num_rows) {
    NML_FERROR(CANNOT_GET_ROW, row, num_rows);
    return NULL;
  }
  nml_mat *r = nml_mat_new(1, num_cols);

  // Good idea: relying on the standard C library to deal with
  //            memory operations
  //memcpy(r->data[0], m->data[row], num_cols * sizeof(*r->data[0]));
  // new form, based in contiguous allocation
  memcpy(r->__data, m->__data + (row * num_cols), num_cols * sizeof(*r->__data));
  return r;
}

// Sets all elements of a matrix to a given value
void nml_mat_all_set(nml_mat *m, double val) {
  const unsigned int num_rows = m->num_rows,
    num_cols = m->num_cols;

  /*int i, j;
  for(i = 0; i < num_rows; i++) {
    for(j = 0; j < num_cols; j++) {
      m->data[i][j] = value;
    }
  }*/

  // Delegate to the well-designed & optimized standard C library, 
  // the task of filling array
  memset( __NML_DATA(m), val, num_rows * num_cols );
}

// Sets all elements of the matrix to given value
int nml_mat_diag_set(nml_mat *m, double value) {
  if (!__NML_IS_SQUARE(m)) {
    NML_FERROR(CANNOT_SET_DIAG, value);
    return R_FAILURE;
  }
  int i;
  const int num_rows = __NML_ROWS(m);
  for(i = 0; i < num_rows; i++) {
    m->data[i][i] = value;

    // assigns:  m[i][j] <-- value
    __NML_ELEM(m, i, i, num_rows) = value;
  }
  return 1;
}

/**
 * Return a new (independent) matrix that is the result of multiplying a given
 * row of the original matrix, by a scalar.
 *
 * @param   m    the original matrix
 * @param   row  the index (zero-based) of the row to multiply
 * @param   num  the scalar
 * @return       the new matrix
 */
nml_mat *nml_mat_row_mult(nml_mat *m, unsigned int row, double num) {
  nml_mat *r = nml_mat_cp(m);
  if (!nml_mat_row_mult_r(r, row, num)) {
    nml_mat_free(r);
    return NULL;
  }
  return r;
}

/**
 * Auxiliary function: Multiply a given row of a matrix, by a scalar. 
 * The matrix will be overwritten.
 *
 * @param   m    the matrix
 * @param   row  index of the row to be multiplied
 * @param   num  the scalar
 * @return       TRUE iff the operation succeeded
 */
bool nml_mat_row_mult_r(nml_mat *m, unsigned int row, double num) {
  const unsigned int num_rows = m->num_rows,
    num_cols = m->num_cols;

  if (row >= num_rows) {
    NML_FERROR(CANNOT_ROW_MULTIPLY, row, num_rows);
    return FALSE;
  }
  int j;
  for(j = 0; j < num_cols; j++) {
    m->data[row][j] *= num;

    // new form, based on contiguous allocation
    __NML_ELEM(m, row, j, num_cols) *= num;
  }
  return TRUE;
}

/**
 * Return a new (independent) matrix that is the result of multiplying a given
 * column of the original matrix, by a scalar.
 *
 * @param   m    the original matrix
 * @param   col  the index (zero-based) of the col to multiply
 * @param   num  the scalar
 * @return       the new matrix
 */
nml_mat *nml_mat_col_mult(nml_mat *m, unsigned int col, double num) {
  nml_mat *r = nml_mat_cp(m);
  if (!nml_mat_col_mult_r(r, col, num)) {
    nml_mat_free(r);
    return NULL;
  }
  return r;
}

/**
 * Auxiliary function: Multiply a given column of a matrix, by a scalar. 
 * The matrix will be overwritten.
 *
 * @param   m    the matrix
 * @param   col  index of the column to be multiplied
 * @param   num  the scalar
 * @return       TRUE iff the operation succeeded
 */
bool nml_mat_col_mult_r(nml_mat *m, unsigned int col, double num) {
  const unsigned int num_rows = __NML_ROWS(m),
    num_cols = __NML_COLS(m);

  if (col >= num_cols) {
    NML_FERROR(CANNOT_COL_MULTIPLY, col, num_cols);
    return FALSE;
  }
  int i;
  for(i = 0; i < num_rows; i++) {
    m->data[i][col] *= num;

    // new form, based on contiguous allocation
    __NML_ELEM(m, i, col, num_cols) *= num;
  }
  return TRUE;
}


/**
 * Elementary row-operation. Adds to the row A[where][:] a multiple of the row A[row][:]
 * This produces an independent matrix.
 */
nml_mat *nml_mat_row_addrow(nml_mat *m, unsigned int where, unsigned int row, double multiplier) {
  nml_mat *r = nml_mat_cp(m);
  if (!nml_mat_row_addrow_r(m, where, row, multiplier)) {
    nml_mat_free(r);
    return NULL;
  }
  return r;
}

/**
 * Elementary row-operation. Adds to the row A[where][:] a multiple of the row A[row][:]
 * The matrix will be overwritten.
 */
int nml_mat_row_addrow_r(nml_mat *m, unsigned int where, unsigned int row, double multiplier) {

  if (where >= m->num_rows || row >= m->num_rows) {
    NML_FERROR(CANNOT_ADD_TO_ROW, multiplier, row, where, m->num_rows);
    return FALSE;
  }
  int j = 0;
  const int num_cols = __NML_COLS(m);
  for(j = 0; j < num_cols; j++) {
    m->data[where][j] += multiplier * m->data[row][j];

    // new form, based on contiguous allocation
    __NML_ELEM(m, where, j, num_cols) += multiplier * __NML_ELEM(m, row, j, num_cols);
  }
  return TRUE;
}

/**
 * Multiply a matrix by a scalar.
 * This produces an independent matrix.
 */
nml_mat *nml_mat_smult(nml_mat *m, double num) {
  nml_mat *r = nml_mat_cp(m);
  nml_mat_smult_r(r, num);
  return r;
}

/**
 * Multiply a matrix by a scalar.
 * Changes are made in-place, the matrix will be overwritten.
 */
int nml_mat_smult_r(nml_mat *m, double num) {
  const unsigned int num_rows = __NML_ROWS(m),
    num_cols = __NML_COLS(m);
  int i, j;
  for(i = 0; i < num_rows; i++) {
    for(j = 0; j < num_cols; j++) {
      m->data[i][j] *= num;

      // new form, based on contiguous allocation
      __NML_ELEM(m, i, j, num_cols) *= num;
    }
  }
  return TRUE;
}

// *****************************************************************************
//
// Modifying the matrix structure
//
// *****************************************************************************

/**
 * Remove column. Produces an independent matrix.
 *
 * @param   m       [description]
 * @param   column  [description]
 * @return          [description]
 */
nml_mat *nml_mat_col_rem(nml_mat *m, unsigned int column) {
  if(column >= m->num_cols) {
    NML_FERROR(CANNOT_REMOVE_COLUMN, column, m->num_cols);
    return NULL;
  }
  nml_mat *r = nml_mat_new(m->num_rows, m->num_cols-1);
  int i, j, k;
  for(i = 0; i < m->num_rows; i++) {
    for(j = 0, k=0; j < m->num_cols; j++) {
      if (column!=j) {
        r->data[i][k++] = m->data[i][j];
      }
    }
  }
  return r;
}

/**
 * Remove row. Produces an independent matrix.
 *
 * @param   m    [description]
 * @param   row  [description]
 * @return       [description]
 */
nml_mat *nml_mat_row_rem(nml_mat *m, unsigned int row) {

  const int num_rows = __NML_ROWS(m), 
    num_cols = __NML_COLS(m);

  if (row >= m->num_rows) {
    NML_FERROR(CANNOT_REMOVE_ROW, row, num_rows);
    return NULL;
  }
  nml_mat *r = nml_mat_new(num_rows-1, num_cols);
  int i, j, k;
  for(i = 0, k = 0; i < num_rows; i++) {
    if (row != i) {
      for(j = 0; j < num_cols; j++) {
        #if 1
        r->data[k][j] = m->data[i][j];
        #endif

        // new form, based on contiguous allocation
        // copying i-th row of m, into k-th row of r
        // r[k][:]  <--  m[i][:]
        memcpy(__NML_DATA(r) + k * num_cols,
          __NML_DATA(m) + i * num_cols,
          num_cols 
        );
      }
      k++;
    }
  }
  return r;
}

/**
 * Swap two rows of a matrix.
 * Produces an independent matrix.
 *
 * @param   m     [description]
 * @param   row1  [description]
 * @param   row2  [description]
 * @return        [description]
 */
nml_mat *nml_mat_row_swap(nml_mat *m, unsigned int row1, unsigned int row2) {
  nml_mat *r = nml_mat_cp(m);
  if (!nml_mat_row_swap_r(r, row1, row2)) {
    nml_mat_free(r);
    return NULL;
  }
  return r;
}

/**
 * This is naive equivalent of std::swap in C++
 */
void nml_swap(double *a, double *b) {
  double tmp = *a;
  *a = *b;
  *b = tmp;
}

/**
 * Swap two rows of a matrix.
 * The matrix will be overwritten.
 *
 * @param   m     [description]
 * @param   row1  [description]
 * @param   row2  [description]
 * @return        [description]
 */
int nml_mat_row_swap_r(nml_mat *m, unsigned int row1, unsigned int row2) {

  if (row1 == row2) { 
    // if row1 == row2, nothing to do
    return TRUE; 
  }

  const int num_rows = __NML_ROWS(m), 
    num_cols = __NML_COLS(m);
  
  if (row1 >= num_rows || row2 >= num_rows) {
    NML_FERROR(CANNOT_SWAP_ROWS, row1, row2, num_rows);
    return FALSE;
  }
  #if 1
  double *tmp = m->data[row2];
  m->data[row2] = m->data[row1];
  m->data[row1] = tmp;
  #endif

  // new form, based on contiguous allocation
  int j;
  for (j = 0; j < num_cols; j++) {
    nml_swap( __NML_DATA(m) + __NML_1D_INDEX(row1, j, num_cols),    // i.e.,  m->__data + (row1*num_cols + j),
                                                                    // but will not incur in overhead as it will be replaced at preprocessing time
      __NML_DATA(m) + __NML_1D_INDEX(row2, j, num_cols)
      );
  }

  return TRUE;
}

/**
 * Swap two columns of a matrix.
 * Produces an independent matrix.
 *
 * @param   m     [description]
 * @param   row1  [description]
 * @param   row2  [description]
 * @return        [description]
 */
nml_mat *nml_mat_col_swap(nml_mat *m, unsigned int col1, unsigned int col2) {
  nml_mat *r = nml_mat_cp(m);
  if (!nml_mat_col_swap_r(r, col1, col2)) {
    nml_mat_free(r);
    return NULL;
  }
  return r;
}

/**
 * Swap two columns of a matrix.
 * Produces an independent matrix.
 *
 * @param   m     [description]
 * @param   row1  [description]
 * @param   row2  [description]
 * @return        [description]
 */
int nml_mat_col_swap_r(nml_mat *m, unsigned int col1, unsigned int col2) {
  
  const int num_rows = __NML_ROWS(m), 
    num_cols = __NML_COLS(m);
  
  if (col1 == col2) { 
    // if col1 == col2, nothing to do
    return TRUE; 
  }

  if (col1 >= num_cols || col2 >= num_rows) {
    NML_FERROR(CANNOT_SWAP_ROWS, col1, col2, num_cols);
    return FALSE;
  }
  double tmp;
  int j;
  #if 1
  for(j = 0; j < num_rows; j++) {
    tmp = m->data[j][col1];
    m->data[j][col1] = m->data[j][col2];
    m->data[j][col2] = tmp;
  }
  #endif

  // new form, based on contiguous allocation
  int i;
  for (i = 0; i < num_rows; i++) {
    nml_swap( __NML_DATA(m) + __NML_1D_INDEX(i, col1, num_cols),    // i.e.,  m->__data + (i*num_cols + col1),
                                                                    // but not incurring in overhead as it will be replaced at preprocessing time
      __NML_DATA(m) + __NML_1D_INDEX(i, col2, num_cols)
      );
  }
  return TRUE;
}

/**
 * Concatenates horizontally a variable number of matrices into one.
 * The concatenation requires the matrices to have the same number of rows, 
 * while the number of columns is allowed to be variable.
 *
 * @param   mnum  number of matrices to be concatenated
 * @param   marr  array of matrices
 * @return        the resulting matrix
 */
nml_mat *nml_mat_cath(unsigned int mnum, nml_mat **marr) {
  if (0 == mnum) {
    return NULL;
  }
  if (1 == mnum) {
    // We just return the one matrix supplied as the first param
    // no need for additional logic
    return nml_mat_cp(marr[0]);
  }
  // We calculate the total number of columns to know how to allocate memory
  // for the resulting matrix
  int i,j,k,offset;
  unsigned int lrow, ncols;
  lrow = marr[0]->num_rows;
  ncols = marr[0]->num_cols;
  for(k = 1; k < mnum; k++) {
    if (NULL == marr[k]) {
      NML_FERROR(INCONSISTENT_ARRAY, k, mnum);
      return NULL;
    }
    if (lrow != marr[k]->num_rows) {
      NML_FERROR(CANNOT_CONCATENATE_H, lrow, marr[k]->num_rows);
      return NULL;
    }
    ncols+=marr[k]->num_cols;
  }
  // At this point we know how the resulting matrix looks like,
  // we allocate memory for it accordingly
  nml_mat *r = nml_mat_new(lrow, ncols);
  #if 1
  for(i = 0; i < r->num_rows; i++) {
    k = 0;
    offset = 0;
    for(j = 0; j < r->num_cols; j++) {
      // If the column index of marr[k] overflows
      if (j-offset == marr[k]->num_cols) {
        offset += marr[k]->num_cols;
        // We jump to the next matrix in the array
        k++;
      }
      r->data[i][j] = marr[k]->data[i][j - offset];
    }
  }
  #endif

  // ______________ THIS NEEDS TO BE TESTED (!!!) _______________
  // new form, based on contiguous allocation
  int pos = 0;    // pos represents the current position into the linear array for r
  for (i = 0; i < lrow; i++) {
    for (k = 0; k < mnum; k++) {    // k-th matrix
      // now, copying the i-th row of the k-th matrix into the linear data for r
      ncols = __NML_COLS(marr[k]);    //number of columns of the current k-th matrix
      memcpy(__NML_DATA(r) + pos, __NML_DATA(marr[k]) + i * ncols, ncols);
      pos += ncols;
    }
  }

  return r;
}

/**
 * Concatenates vertically a variable number of matrices into one.
 * The concatenation requires the matrices to have the same number of columns, 
 * while the number of rows is allowed to be variable.
 *
 * @param   mnum  number of matrices to be concatenated
 * @param   marr  array of matrices
 * @return        the resulting matrix
 */
nml_mat *nml_mat_catv(unsigned int mnum, nml_mat **marr) {
  if (0 == mnum) {
    return NULL;
  }
  if (1 == mnum) {
    return nml_mat_cp(marr[0]);
  }
  // We check to see if the matrices have the same number of columns
  int lcol, i, j, k, offset;
  unsigned int numrows;
  nml_mat *r;
  lcol = marr[0]->num_cols;
  numrows = 0;
  for(i = 0; i < mnum; i++) {
    if (NULL==marr[i]) {
      NML_FERROR(INCONSISTENT_ARRAY, i, mnum);
      return NULL;
    }
    if (lcol != marr[i]->num_cols) {
      NML_FERROR(CANNOT_CONCATENATE_V,lcol,marr[i]->num_cols);
      return NULL;
    }
    // In the same time we calculate the resulting matrix number of rows
    numrows+=marr[i]->num_rows;
  }
  // At this point we know the dimensions of the resulting Matrix
  r = nml_mat_new(numrows, lcol);
  // We start copying the values one by one
  #if 1
  for(j = 0; j < r->num_cols; j++) {
    offset = 0;
    k = 0;
    for(i = 0; i < r->num_rows; i++) {
      if (i - offset == marr[k]->num_rows) {
        offset += marr[k]->num_rows;
        k++;
      }
      r->data[i][j] = marr[k]->data[i-offset][j];
    }
  }
  #endif

  // ______________ THIS NEEDS TO BE TESTED (!!!) _______________
  // new form, based on contiguous allocation
  int pos = 0;    // pos represents the current position into the linear array for r
  for (k = 0; k < mnum; k++) {    // k-th matrix
    for (i = 0; i < __NML_ROWS(marr[k]); i++) {
      // now, copying the i-th row of the k-th matrix into the linear data for r
      memcpy(__NML_DATA(r) + pos, __NML_DATA(marr[k]) + i * lcol, lcol);
      pos += lcol;
    }
  }

  nml_mat_print(r);
  return r;
}

// *****************************************************************************
//
// Matrix Operations
//
// *****************************************************************************
//
// Matrix Operations
//

/**
 * Adding matrices.
 *
 * @param   m1  matrix 1
 * @param   m2  matrix 2
 * @return      a pointer for the new matrix (m1 + m2)
 */
nml_mat *nml_mat_add(nml_mat *m1, nml_mat *m2) {
  nml_mat *r = nml_mat_cp(m1);
  if (!nml_mat_add_r(r, m2)) {
    nml_mat_free(r);
    return NULL;
  }
  return r;
}

/**
 * Adding matrices.
 *
 * @param   m1  matrix 1
 * @param   m2  matrix 2
 * @return      After operation, m1 will contain (m1 + m2).
 */
int nml_mat_add_r(nml_mat *m1, nml_mat *m2) {
  if (!nml_mat_eqdim(m1, m2)) {
    NML_ERROR(CANNOT_ADD);
    return FALSE;
  }
  int i, j;
  #if 1
  for(i = 0; i < m1->num_rows; i++) {
    for(j = 0; j < m1->num_cols; j++) {
      m1->data[i][j] += m2->data[i][j];
    }
  }
  #endif

  // new form, based on contiguous allocation
  int k;
  double *p1 = __NML_DATA(m1);
  double *p2 = __NML_DATA(m2);
  for (k = 0; k < __NML_ROWS(m1) * __NML_COLS(m1); k++) {
    *p1 += *p2;
  }
  return TRUE;
}

/**
 * Substracting matrices.
 *
 * @param   m1  matrix 1
 * @param   m2  matrix 2
 * @return      a pointer for the new matrix (m1 - m2)
 */
nml_mat *nml_mat_sub(nml_mat *m1, nml_mat *m2) {
  nml_mat *r = nml_mat_cp(m2);
  if (!nml_mat_sub_r(r, m2)) {
    nml_mat_free(r);
    return NULL;
  }
  return r;
}

/**
 * Substracting matrices.
 *
 * @param   m1  matrix 1
 * @param   m2  matrix 2
 * @return      After operation, m1 will contain (m1 - m2).
 */
int nml_mat_sub_r(nml_mat *m1, nml_mat *m2) {
  if (!nml_mat_eqdim(m1, m2)) {
    NML_ERROR(CANNOT_SUBSTRACT);
    return 0;
  }
  #if 1
  int i, j;
  for(i = 0; i < m1->num_rows; i++) {
    for(j = 0; j < m1->num_cols; j++) {
      m1->data[i][j] -= m2->data[i][j];
    }
  }
  #endif

  // new form, based on contiguous allocation
  int k;
  double *p1 = __NML_DATA(m1);
  double *p2 = __NML_DATA(m2);
  for (k = 0; k < __NML_ROWS(m1) * __NML_COLS(m1); k++) {
    *p1 -= *p2;
  }
  return 1;
}

/**
 * Matrix product.
 *
 * @param   m1  matrix 1
 * @param   m2  matrix 2
 * @return      a pointer for the new matrix (m1 * m2)
 */
nml_mat *nml_mat_dot(nml_mat *m1, nml_mat *m2) {
  if (!(m1->num_cols == m2->num_rows)) {
    NML_ERROR(CANNOT_MULTIPLY);
    return NULL;
  }
  int i, j, k;
  nml_mat *r = nml_mat_new(m1->num_rows, m2->num_cols);
  int  r_ncols = __NML_COLS(r);
  int m1_ncols = __NML_COLS(m1);
  int m2_ncols = __NML_COLS(m2);
  for(i = 0; i < r->num_rows; i++) {
    for(j = 0; j < r->num_cols; j++) {
      for(k = 0; k < m1->num_cols; k++) {
        r->data[i][j] += m1->data[i][k] * m2->data[k][j];

        // new form, based on contiguous allocation
        // r->data[i][j] += m1->data[i][k] * m2->data[k][j];
        __NML_ELEM(r, i, j, r_ncols) += __NML_ELEM(m1, i, k, m1_ncols) * 
          __NML_ELEM(m2, k, j, m2_ncols);
      }
    }
  }
  return r;
}

/**
 * Transpose matrix
 *
 * @param   m  the matrix
 * @return     pointer to the created transpose of m
 */
nml_mat *nml_mat_transp(nml_mat *m) {
  int i, j;
  const int num_rows = __NML_ROWS(m), 
    num_cols = __NML_COLS(m);

  nml_mat *r = nml_mat_new(num_cols, num_rows);
  for(i = 0; i < r->num_rows; i++) {
    for(j = 0; j < r->num_cols; j++) {
      r->data[i][j] = m->data[j][i];

      // new form, based on contiguous allocation
      // note: r->num_cols is m->num_rows
      __NML_ELEM(r, i, j, num_rows) = __NML_ELEM(m, j, j, num_cols);
    }
  }
  return r;
}

/**
 * Trace of the matrix
 *
 * @param   m  the matrix
 * @return     trace of m
 */
double nml_mat_trace(nml_mat* m) {
  const unsigned int num_cols = __NML_COLS(m);

  if (!m->is_square) {
    NML_ERROR(CANNOT_TRACE);
  }
  int i;
  double trace = 0.0;
  for(i = 0; i < m->num_rows; i++) {
    trace += m->data[i][i];

    // new form, based on contiguous allocation
    trace += __NML_ELEM(m, i, i, num_cols);
  }
  return trace;
}

// *****************************************************************************
//
// Row Echelon
//
// *****************************************************************************

/**
 * Find the first non-zero element in the column `col`, and on or under the row `row`.
 * This is used to determine the pivot in Gauss elimination.
 * If no pivot is found, returns -1
 */
int _nml_mat_pivotidx(nml_mat *m, unsigned int col, unsigned int row) {
  // No validations are made, this is an API Method
  int i;
  const unsigned int num_rows = __NML_ROWS(m),
    num_cols = __NML_COLS(m);
  for (i = row; i < num_rows; i++) {
    //if (fabs(m->data[i][col]) > NML_MIN_COEF) {
    if (__NML_ELEM(m, i, col, num_cols) > NML_MIN_COEF) {
      return i;
    }
  }
  return -1;
}

/**
 * Find the element with the max absolute value in the column `col`, 
 * and on or under the row `row`.
 * This is used to determine the pivot in Gauss elimination.
 * If no pivot is found, returns -1
 */
int _nml_mat_pivotmaxidx(nml_mat *m, unsigned int col, unsigned int row) {
  int i, max_i;
  double pivot;     //micol;     // `micol`?   maybe you mean `pivot`?
  double max = fabs(m->data[row][col]);
  const int num_rows = __NML_ROWS(m), 
    num_cols = __NML_COLS(m);

  max_i = row;
  for(i = row + 1; i < num_rows; i++) {    // this should be for i > row, only
    //pivot = fabs(m->data[i][col]);
    pivot = fabs(__NML_ELEM(m, i, col, num_cols));    // i.e., fabs( m[i, col] )
    if (pivot > max) {
      max = pivot;
      max_i = i;
    }
  }
  return (max < NML_MIN_COEF) ? -1 : max_i;
}

/**
 * Retrieves the matrix in Row Echelon form using Gauss Elimination
 *
 * @param   m  the argument matrix
 * @return     the matrix reduced to its Row Echelon form
 */
nml_mat *nml_mat_ref(nml_mat *m) {

  const int num_cols = __NML_COLS(m);
  int i, j, k, pivot;
  double coef;

  i = j = 0;
  nml_mat *r = nml_mat_cp(m);
  while (j < num_cols && i < num_cols) {
    // Find the pivot - the first non-zero entry in the first column of the matrix
    //pivot = _nml_mat_pivotidx(r, j, i);
    pivot = _nml_mat_pivotmaxidx(r, j, i);
    if (pivot<0) {
      // All elements on the column are zeros
      // We move to the next column without doing anything
      j++;
      continue;
    }
    // We interchange rows moving the pivot to the first row that doesn't have
    // already a pivot in place
    if (pivot != i) {
      nml_mat_row_swap_r(r, i, pivot);
    }
    // Multiply each element in the pivot row by the inverse of the pivot
    //nml_mat_row_mult_r(r, i, 1/r->data[i][j]);
    nml_mat_row_mult_r(r, i, 1.0 / __NML_ELEM(r, i, j, num_cols));
    
    // We add multiplies of the pivot so every element on the column equals 0
    for(k = i+1; k < r->num_rows; k++) {
      coef = __NML_ELEM(r, k, j, num_cols);
      if (fabs(coef) > NML_MIN_COEF) {
        nml_mat_row_addrow_r(r, k, i, -coef);
      } 
    }
    i++;
    j++;
  }
  return r;
}

/** 
 * Retrieves the matrix in Reduced Row Echelon using Gauss-Jordan Elimination 
 */
nml_mat *nml_mat_rref(nml_mat *m) {
  
  const int num_rows = __NML_ROWS(m),
    num_cols = __NML_COLS(m);
  int i, j, k, pivot;

  nml_mat* r = nml_mat_cp(m);
  i = j = 0;
  while (j < num_cols && i < num_rows) {
    // We find the pivot, the maximum row id (fabs) in the column
    pivot = _nml_mat_pivotmaxidx(r, j, i);
    if (pivot < 0) {
      // No pivot, we change columns
      j++;
      continue;
    }
    // We interchange rows to out the pivot row into the 
    // desired position
    if (pivot != i) {
      nml_mat_row_swap_r(r, i, pivot);
    }
    // We create 1 in the pivot position
    nml_mat_row_mult_r(r, i, 1.0 / __NML_ELEM(r, i, j, num_cols));
     // We put zeros on the colum with the pivot
    for (k = 0; k < num_rows; k++) {
      if (k != i) {
        nml_mat_row_addrow_r(r, k, i, - __NML_ELEM(r, k, j, num_cols) );
      }
    }
    i++;
    j++;
  }
  return r;
}

// *****************************************************************************
//
// LUP Decomposition
//
// *****************************************************************************

// Finds the maxid on the column (starting from k -> num_rows)
// This method is used for pivoting in LUP decomposition
int _nml_mat_absmaxr(nml_mat *m, unsigned int k) {
  // Find max id on the column;
  int i;
  double max = m->data[k][k];
  int maxIdx = k;
  for(i = k+1; i < m->num_rows; i++) {
    if (fabs(m->data[i][k]) > max) {
      max = fabs(m->data[i][k]);
      maxIdx = i;
    }
  }
  return maxIdx;
}

/**
 * Yoel.- [2021.04.28]
 * I believe, that `_nml_mat_absmaxr(m,k)` will produce the same than
 * `_nml_mat_pivotmaxidx(m,k,k)`
 */

/**
 * Allocates memory for a new nml_mat_lup structure
 *
 * @param   L                 lower triangular matrix
 * @param   U                 upper triangular matrix
 * @param   P                 permutation matrix
 * @param   num_permutations  number of permutations
 * @return                    pointer to LUP structure allocated
 */
nml_mat_lup *nml_mat_lup_new(nml_mat *L, nml_mat *U, nml_mat *P, unsigned int num_permutations) {
  nml_mat_lup *r = malloc(sizeof(*r));
  NP_CHECK(r);
  r->L = L;
  r->U = U;
  r->P = P;
  r->num_permutations = num_permutations;
  return r;
}

/**
 * Free a LUP structure
 */
void nml_mat_lup_free(nml_mat_lup* lu) {
  nml_mat_free(lu->P);
  nml_mat_free(lu->L);
  nml_mat_free(lu->U);
  free(lu);
}

/**
 * Print a LUP structure
 */
void nml_mat_lup_print(nml_mat_lup *lu) {
  nml_mat_print(lu->L);
  nml_mat_print(lu->U);
  nml_mat_print(lu->P);
}

/**
 * Print a LUP structure, with format
 */
void nml_mat_lup_printf(nml_mat_lup *lu, const char *fmt) {
  nml_mat_printf(lu->L, fmt);
  nml_mat_printf(lu->U, fmt);
  nml_mat_printf(lu->P, fmt);
}

/**
 * LU(P) factorization of a matrix.
 *
 * @param   m  the matrix
 * @return     pointer to a LUP structure.
 */
nml_mat_lup *nml_mat_lup_solve(nml_mat *m) {

  const int num_rows = __NML_ROWS(m),
    num_cols = __NML_COLS(m);
  int j,i, pivot;
  unsigned int num_permutations = 0;
  double lambda;
  
  if ( __NML_IS_SQUARE(m) ) {
    NML_FERROR(CANNOT_LU_MATRIX_SQUARE, num_rows, num_cols);
    return NULL;
  }
  nml_mat *L = nml_mat_new(num_rows, num_rows);
  nml_mat *U = nml_mat_cp(m);
  nml_mat *P = nml_mat_eye(num_rows);

  for(j = 0; j < num_cols; j++) {
    // Retrieves the row with the biggest element for column (j)
    pivot = _nml_mat_absmaxr(U, j);
    if (fabs(U->data[pivot][j]) < NML_MIN_COEF) {
      NML_ERROR(CANNOT_LU_MATRIX_DEGENERATE);
      return NULL;
    }
    if (pivot != j) {
      // Pivots LU and P accordingly to the rule
      nml_mat_row_swap_r(U, j, pivot);
      nml_mat_row_swap_r(L, j, pivot);
      nml_mat_row_swap_r(P, j, pivot);
      // Keep the number of permutations to easily calculate the
      // determinant sign afterwards
      num_permutations++;
    }
    for(i = j+1; i < num_rows; i++) {
      //lambda = U->data[i][j] / U->data[j][j];    // lambda: multiplier
      //
      // new form
      // lambda: this is the multiplier for making U[i,j] to be zero
      lambda = __NML_ELEM(U, i, j, num_cols) / __NML_ELEM(U, j, j, num_cols);

      // Building the U upper rows
      nml_mat_row_addrow_r(U, i, j, -lambda);
      // Store the multiplier in L
      //L->data[i][j] = lambda;
      //
      //new form
      __NML_ELEM(L, i, j, num_rows) = lambda;    
      //                  ~~~~~~~~
      //                       ^_________________ L is (m x m), if  is (m x n)
      //     ^            
      //     |___________________________________ replaces L->data[i][j]
    }
  }
  nml_mat_diag_set(L, 1.0);

  return nml_mat_lup_new(L, U, P, num_permutations);
}

// After the LU(P) factorisation the determinant can be easily calculated
// by multiplying the main diagonal of matrix U with the sign.
// the sign is -1 if the number of permutations is odd
// the sign is +1 if the number of permutations is even
double nml_mat_det(nml_mat_lup* lup) {
  const int num_rows = __NML_ROWS(lup->U),
    num_cols = __NML_COLS(lup->U);
  int k;

  int sign = (lup->num_permutations%2==0) ? 1 : -1;
  nml_mat *U = lup->U;
  double product = 1.0;
  for(k = 0; k < num_rows; k++) {
    //product *= U->data[k][k];
    //
    // new form
    product *= __NML_ELEM(U, k, k, num_cols);
    //         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                      ^
    //                      |_______________ U[k][k]
  }
  return product * sign;
}

// Returns LU matrix from a LUP structure
nml_mat *nml_mat_lu_get(nml_mat_lup* lup) {
  nml_mat *r = nml_mat_cp(lup->U);
  // Copy L (without first diagonal in result)
  int i, j;
  nml_mat *L = lup->L;
  for(i = 1; i < L->num_rows; i++) {
    for(j = 0; j < i; j++) {
      //r->data[i][j] = L->data[i][j];
      //
      // new form
      __NML_ELEM(r, i, j, L->num_cols) = __NML_ELEM(L, i, j, L->num_cols);
    }
  }
  return r;
}

// *****************************************************************************
//
// Solving linear systems of equations
//
// *****************************************************************************

// Forward substitution algorithm.
// Solves the linear system L*x = b
//
// L is lower triangular matrix of size NxN
// b is column matrix of size Nx1
// x is the solution column matrix of size Nx1
//
// Note: In case L is not a lower triangular matrix, the algorithm will try to
// select only the lower triangular part of the matrix L and solve the system
// with it.
//
// Note: In case any of the diagonal elements (L[i][i]) are 0 the system cannot
// be solved
//
// Note: This function is usually used with an L matrix from a LU decomposition
nml_mat *nml_ls_solvefwd(nml_mat *L, nml_mat *b) {
  const int num_rows = __NML_ROWS(L),
    num_cols = __NML_COLS(L);

  // assuring L is square
  if (num_rows != num_cols) {
    NML_FERROR(MATRIX_NOT_SQUARE, num_rows, num_cols);
  }

  nml_mat* x = nml_mat_new(num_cols, 1);
  int i, j;
  double tmp, diag;
  for(i = 0; i < num_cols; i++) {

    /* 
     * x[i] = ( b[i] -  sum  { L[i][j] * x[j] } ) / L[i][i]
     *                 j < i
     */
    #if NML_V1    // version 1
    tmp = b->data[i][0];
    for(j = 0; j < i ; j++) {
      tmp -= L->data[i][j] * x->data[j][0];
    }
    x->data[i][0] = tmp / L->data[i][i];
    #endif

    // now, traslating into version 2, with the new form
    #if NML_V2
    tmp = __NML_ELEM2(b, i, 0);    // b[i];
    for(j = 0; j < i ; j++) {
      tmp -= __NML_ELEM2(L, i, j) * __NML_ELEM2(x, j, 0);    //  -  sum  { L[i][j] * x[j] }
                                                             //   (j < i)
    }

    // assuring diagonal element is not zero
    diag = __NML_ELEM2(L, i, i);
    if (fabs(diag) < NML_MIN_COEF) {
      NML_ERROR(INCONSISTENT_SYSTEM);
      return NULL;
    }
    __NML_ELEM2(x, i, 0) = tmp / diag;    // x[i] = tmp / L[i][i]
    #endif    
  }
  return x;
}


// Back substitution algorithm.
// Solves the linear system U*x = b
//
// U is an upper triangular matrix of size NxN
// b is a column matrix of size Nx1
// x is the solution column matrix of size Nx1
//
// Note in case U is not an upper triangular matrix, the algorithm will try to
// select only the upper triangular part of the matrix U and solve the system
// with it
//
// Note: In case any of the diagonal elements (U[i][i]) are 0 the system cannot
// be solved
nml_mat *nml_ls_solvebck(nml_mat *U, nml_mat *b) {
  const int num_rows = __NML_ROWS(U),
    num_cols = __NML_COLS(U);

  // assuring L is square
  if (num_rows != num_cols) {
    NML_FERROR(MATRIX_NOT_SQUARE, num_rows, num_cols);
  }

  nml_mat *x = nml_mat_new(U->num_cols, 1);
  int i, j;
  double tmp, diag;
  for (i = num_cols - 1; i >= 0; i--) {
    /* 
     * x[i] = ( b[i] -  sum  { U[i][j] * x[j] } ) / U[i][i]
     *                 j > i
     */
    #if NML_V1    // version 1
    tmp = b->data[i][0];
    for(j = i; j < U->num_cols; j++) {
      tmp -= U->data[i][j] * x->data[j][0];
    }
    x->data[i][0] = tmp / U->data[i][i];
    #endif

    // now, traslating into version 2, with the new form
    #if NML_V2
    tmp = __NML_ELEM2(b, i, 0);    // b[i];
    for(j = 0; j < i ; j++) {
      tmp -= __NML_ELEM2(U, i, j) * __NML_ELEM2(x, j, 0);    //  -  sum  { U[i][j] * x[j] }
                                                             //   (j < i)
    }

    // assuring diagonal element is not zero
    diag = __NML_ELEM2(U, i, i);    // U[i][i]
    if (fabs(diag) < NML_MIN_COEF) {
      NML_ERROR(INCONSISTENT_SYSTEM);
      return NULL;
    }
    __NML_ELEM2(x, i, 0) = tmp / diag;
    #endif    
  }
  return x;
}

/**
 * A[n][n] is a square matrix
 * m contains matrices L, U, P for A[n][n] so that P*A = L*U
 * 
 * The linear system is:
 * 
 *   A*x=b  =>  P*A*x = P*b  =>  L*U*x = P*b 
 *   
 * (where b is a matrix[n][1], and x is a matrix[n][1])
 * if y = U*x , we solve two systems:
 *    L * y = P b (forward substition)
 *    U * x = y (backward substition)
 *    
 * We obtain and return x.
 */
nml_mat *nml_ls_solve(nml_mat_lup *lu, nml_mat* b) {
  if (lu->U->num_rows != b->num_rows || b->num_cols != 1) {
    NML_FERROR(CANNOT_SOLVE_LIN_SYS_INVALID_B,
      b->num_rows,
      b->num_cols,
      lu->U->num_rows,
      1);
      return NULL;
  }
  nml_mat *Pb = nml_mat_dot(lu->P, b);

  // We solve L*y = P*b using forward substition
  nml_mat *y = nml_ls_solvefwd(lu->L, Pb);

  // We solve U*x=y
  nml_mat *x = nml_ls_solvebck(lu->U, y);

  nml_mat_free(y);
  nml_mat_free(Pb);
  return x;
}

/**
 * Calculates the inverse of a matrix. To find the inverse B of A,
 * we solve the linear system AB = I, where I is the identity matrix
 * (of same size of A). The system is solved by LUP decomposition.
 *
 * @param   lup  the LU decomposition of A.
 * @return       the inverse of A.
 */
nml_mat *nml_mat_inv(nml_mat_lup *lup) {
  const unsigned int num_cols = lup->L->num_cols;
  nml_mat *r = nml_mat_sqr(num_cols);
  nml_mat *I = nml_mat_eye(lup->U->num_rows);
  nml_mat *invx;
  nml_mat *Ix;
  int i,j;
  for(j = 0; j < num_cols; j++) {
    Ix = nml_mat_col_get(I, j);
    invx = nml_ls_solve(lup, Ix);

    // the j-th column of the A^(-1) is the column vector result
    for(i = 0; i < invx->num_rows; i++) {
      #if NML_V1    // version 1
      r->data[i][j] = invx->data[i][0];
      #endif

      #if NML_V2    // version 2
      __NML_ELEM2(r, i, j) = __NML_ELEM(invx, i, 0, 1);
      #endif

    }
    nml_mat_free(invx);
    nml_mat_free(Ix);
  }
  nml_mat_free(I);
  return r;
}

// *****************************************************************************
//
// QR Decomposition
//
// *****************************************************************************

// Useful for QR decomposition
// Represents the (dot) product of two vectors:
// vector1 = m1col column from m1
// vector2 = m2col column from m2
double nml_vect_dot(nml_mat *m1, unsigned int m1col, nml_mat *m2, unsigned m2col) {
  if (m1->num_rows!=m2->num_rows) {
    NML_FERROR(CANNOT_VECT_DOT_DIMENSIONS, m1->num_rows, m2->num_rows);
  }
  if (m1col >= m1->num_cols) {
    NML_FERROR(CANNOT_GET_COLUMN, m1col, m1->num_cols);
  }
  if (m2col >= m2->num_cols) {
    NML_FERROR(CANNOT_GET_COLUMN, m2col, m2->num_cols);
  }
  int i;
  double dot = 0.0;
  const unsigned num_cols1 = m1->num_cols;
  const unsigned num_cols2 = m2->num_cols;
  for(i = 0; i < m1->num_rows; i++) {
    #if NML_V1    // NML version 1
    dot += m1->data[i][m1col] * m2->data[i][m2col];
    #endif

    #if NML_V2    // NML version 2
    dot += __NML_ELEM(m1, i, m1col, num_cols1) * __NML_ELEM(m2, i, m2col, num_cols2);
    #endif
  }
  return dot;
}

/**
 * Calculates the l2 norm for a colum in the matrix
 *
 * @param   m    the matrix
 * @param   col  index of column for which to calculate the l2 norm
 * @return       l2 norm of column `col`
 */
double nml_mat_col_l2norm(nml_mat *m, unsigned int col) {
  const unsigned num_rows = __NML_ROWS(m), 
    num_cols = __NML_COLS(m);
  double square_sum = 0.0;
  int i;

  if(col >= num_cols) {
    NML_FERROR(CANNOT_COLUMN_L2NORM, col, num_cols);
  }
  for(i = 0; i < num_rows; i++) {
    #if NML_V1    // NML version 1
    square_sum += (m->data[i][col]*m->data[i][col]);
    #endif

    /**
     * l2norm(col) = (    sum  { m[i][col] ^ 2 }  ) ^ (1/2)
     *                1 <= i <= n
     */
    #if NML_V2    // NML version 2
    square_sum += (__NML_ELEM2(m, i, col) * __NML_ELEM2(m, i, col));
    /*             ~~~~~~~~~~~~~~~~~~~~~~
     *                                  ^ this definition takes the constant `num_cols`
     *                                  | as the 2nd dimension of m
     */
    #endif
  }
  return sqrt(square_sum);
}

/**
 * Calculates the l2norm for each column
 * Keeps results into 1 row matrix, where the j-th element
 * corresponds to the l2norm for the j-th column of m
 *
 * @param   m  the matrix
 * @return     row matrix, containing the l2norm at colum-wise.
 */
nml_mat *nml_mat_l2norm(nml_mat *m) {
  const unsigned num_rows = __NML_ROWS(m), 
    num_cols = __NML_COLS(m);
  int i, j;
  double square_sum;
  
  nml_mat *r = nml_mat_new(1, num_cols);
  for(j = 0; j < num_cols; j++) {
    square_sum = 0.0;
    #if NML_V1    // NML version 1
    for(i = 0; i < num_rows; i++) {
      square_sum += m->data[i][j] * m->data[i][j];
    }
    r->data[0][j] = sqrt(square_sum);
    #endif

    /**
     * l2norm[j] = (    sum  { m[i][j] ^ 2 }  ) ^ (1/2)
     *              1 <= i <= n
     */
    #if NML_V2    // NML version 2
    for(i = 0; i < num_rows; i++) {
      square_sum += m->data[i][j] * m->data[i][j];
      square_sum += (__NML_ELEM2(m, i, j) * __NML_ELEM2(m, i, j));
      /*             ~~~~~~~~~~~~~~~~~~~~
       *                                ^ this definition takes the constant `num_cols`
       *                                | as the 2nd dimension of m
       */
    }
    __NML_ELEM2(r, 0, j) = sqrt(square_sum);
    #endif
  }
  return r;
}

/**
 * Normalizes a matrix, that is, divides each column by its own l2norm.
 * This way, each column will have an unitary square-norm.
 *
 * @param   m  the matrix
 * @return     the normalized matrix (using the same storage space than 
 *             its predecesor m)
 */
nml_mat *nml_mat_normalize(nml_mat *m) {
  nml_mat *r = nml_mat_cp(m);
  if (!nml_mat_normalize_r(r)) {
    nml_mat_free(r);
    return NULL;
  }
  return r;
}

/**
 * Auxiliary function to normalize matrix m.
 */
int nml_mat_normalize_r(nml_mat *m) {
  const unsigned num_rows = __NML_ROWS(m), 
    num_cols = __NML_COLS(m);
  int j;

  nml_mat *l2norms = nml_mat_l2norm(m);
  for(j = 0; j < m->num_cols; j++) {
    #if NML_V1    // NML version 1
    if (l2norms->data[0][j] < NML_MIN_COEF) {
      NML_FERROR(VECTOR_J_DEGENERATE, j);
      nml_mat_free(l2norms);
      return 0;
    }
    nml_mat_col_mult_r(m, j, 1/l2norms->data[0][j]);
    #endif

    #if NML_V2    // NML version 2
    if (__NML_ELEM2(l2norms, 0, j) < NML_MIN_COEF) {    // l2norms[0][j]
      NML_FERROR(VECTOR_J_DEGENERATE, j);
      nml_mat_free(l2norms);
      return 0;
    }
    nml_mat_col_mult_r(m, j, 1 / __NML_ELEM2(l2norms, 0, j));    // l2norms[0][j]
    #endif
  }
  nml_mat_free(l2norms);
  return 1;
}

/**
 * Allocates a new QR structure.
 *
 * @return  new QR structure
 */
nml_mat_qr *nml_mat_qr_new() {
  nml_mat_qr *qr = malloc(sizeof(*qr));
  NP_CHECK(qr);
  return qr;
}

/**
 * Deallocates a QR structure.
 *
 * @param  qr  an existing, previously allocated, QR structure.
 */
void nml_mat_qr_free(nml_mat_qr *qr) {
  nml_mat_free(qr->Q);
  nml_mat_free(qr->R);
  free(qr);
}

/**
 * QR factorization. The QR structure will contain both, the ortoghonal matrix Q,
 * and the upper triangular matrix R, so that M = QR.
 *
 * @param   m  the matrix
 * @return     QR structure containing the matrices Q and R.
 */
nml_mat_qr *nml_mat_qr_solve(nml_mat *m) {
  const unsigned num_rows = __NML_ROWS(m), 
    num_cols = __NML_COLS(m);
  int j, k;
  double l2norm;
  double rkj;

  // creating QR structure
  nml_mat_qr *qr = nml_mat_qr_new();
  nml_mat *Q = nml_mat_cp(m);
  nml_mat *R = nml_mat_new(num_rows, num_cols);
  nml_mat *aj;
  nml_mat *qk;
  // Gram-Schmidt process
  for(j = 0; j < num_cols; j++) {    
    rkj = 0.0;
    aj = nml_mat_col_get(m, j);
    #if NML_V1    // NML version 1
    for(k = 0; k < j; k++) {
       rkj = nml_vect_dot(m, j, Q, k);
       R->data[k][j] = rkj;
       qk = nml_mat_col_get(Q, k);
       nml_mat_col_mult_r(qk, 0, rkj);
       nml_mat_sub_r(aj, qk);
       nml_mat_free(qk);
    }
    for(k = 0; k < num_rows; k++) {
      Q->data[k][j] = aj->data[k][0];
    }
    l2norm = nml_mat_col_l2norm(Q, j);
    nml_mat_col_mult_r(Q, j, 1/l2norm);
    R->data[j][j] = l2norm;
    #endif

    #if NML_V2    // NML version 2
    for(k = 0; k < j; k++) {
      // dot product of j-th column of m, and k-th column of Q
      rkj = nml_vect_dot(m, j, Q, k);
      __NML_ELEM2(R, k, j) = rkj;    // R[k][j] = rkj
      /* ~~~~~~~~~~~~~~~~~
       *                 ^ this definition takes constant `num_cols` as the 2nd dimension of R
       */
      
      /* M_j = M_j -  sum { dot(M_j, Q_k) }
       *             k < j  
       */
      qk = nml_mat_col_get(Q, k);
      nml_mat_col_mult_r(qk, 0, rkj);
      nml_mat_sub_r(aj, qk);
      nml_mat_free(qk);
    }
    for(k = 0; k < num_rows; k++) {
      Q->data[k][j] = aj->data[k][0];
    }
    l2norm = nml_mat_col_l2norm(Q, j);
    nml_mat_col_mult_r(Q, j, 1/l2norm);
    R->data[j][j] = l2norm;
    #endif    
    nml_mat_free(aj);
  }
  qr->Q = Q;
  qr->R = R;
  return qr;
}

#ifdef __COLS
#undef __COLS 
#endif