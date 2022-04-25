/** 
 * NML_DEF_H_
 *
 * Additional definitions.
 * 
 * Added by Yoel.- 2022.04.25
 */

#ifndef NML_DEF_H_
#  define NML_DEF_H_

/* Return status'es */
#ifndef R_SUCCESS
      #define R_SUCCESS 0
#endif
#ifndef R_FAILURE
      #define R_FAILURE 1
#endif

/* boolean types */
#ifndef bool
enum __enum_bool {FALSE=0, TRUE};
typedef enum __enum_bool bool_t;
typedef bool_t bool;
#endif

#ifndef True    // aliases for TRUE, FALSE
      #define True TRUE
#endif
#ifndef true
      #define true TRUE
#endif
#ifndef False
      #define False FALSE
#endif
#ifndef false
      #define false FALSE
#endif

#define _max(a,b) ((a) >= (b) ? a : b)

/**
 * Yoel.- 2022-04-25
 *
 * This is to simplify the nomenclature for accessing the structures' data members
 * 
 * Inspired in the gnu stdc name convention:
 * https://github.com/gcc-mirror/gcc/blob/d9375e490072d1aae73a93949aa158fcd2a27018/libstdc%2B%2B-v3/include/bits/stl_algo.h#L3858
 */
#define __COLS(m) ((m)->num_cols)
#define __ROWS(m) ((m)->num_rows)
#define __DATA(m) ((m)->__data)

/**
 * this MACRO (not function) converts 2D indexes into a 1D contiguous index, by means of
 *
 * [i,j] -> i * n + j
 *
 * @param   i  row index
 * @param   j  column index
 * @param   n  2nd dimension: number of columns
 */
#define __1D_INDEX(i, j, n) ((i)*(n) + (j))

/**
 * It is equivalent to matrix->data[i][j], in 2D indexing
 */
#define __ELEM(matrix, i, j) ((matrix)->__data[ __1D_INDEX(i, j, __COLS(matrix)) ])

#endif