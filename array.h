#ifndef __ARRAY_H__
#define __ARRAY_H__

#include "common.h"
#include <stdint.h>
#include <stdlib.h>

#ifdef __ARRAY_MAIN__
#define EXT_ARRAY
#else
#define EXT_ARRAY extern
#endif

typedef enum { 
    undefArray, int8Array, int16Array, int32Array, int64Array, floatArray, doubleArray 
} ArrayType;

EXT_ARRAY char    const * const arrayTypeStrings[];      // string names of each type
EXT_ARRAY size_t  const arrayTypeElementSizes[]; // size of one element in bytes 

typedef struct ArrayStruct *Array;

Array Array_new(ArrayType type);                 // initialize a new array 
                                                 // instead of Array x; always use Array x = Array_new(<type>)
void Array_setType(Array a, ArrayType type);     // set the type later if not known at variable creation.
ArrayType Array_getType(Array a);         

/* Array allocation */
bool Array_allocate_dims(Array a, int ndims, int *dims); // malloc dims[0]*...*dims[ndims]
bool Array_allocate(Array a, int n);                     // malloc n*sizeof(type), ndims <-1, dimsize[1] <- n
bool Array_setDimensions(Array a, int ndims, int *dims); // use only after simple allocate
                                                         // dims[0]*...*dims[ndims] must = Array_length(a)
void Array_free(Array a);

/* Array information */
long Array_length(Array a);                      // = n
int Array_getDimensions(Array a);                // return ndims, 
int Array_getDimensionSize(Array a, int dimidx); // return dimsize[dimidx] or -1 if idx is out of scope
int Array_getAllDimensions(Array a, int *dims);   // return ndims, and fill dims[] with the values. no alloc inside!
                                                 // if dims is NULL, only the size is returned
void Array_calcminmax(Array a);                  // calculate the min and max. Use before asking for min or max
double Array_min(Array a);
double Array_max(Array a);
int Array_errno(Array a);                        // 0 if last op was okay
void Array_seterrno(Array a, int err);           // readers can set error state
void Array_strvalue(Array a, long idx, char* str); // return in str the string value of a[idx]
void Array_strvalue_padded(Array a, long idx, char* str, int width); // pad with zeros to given width
double Array_doublevalue(Array a, long idx);      // return value of a[idx] as double

/** Direct access to the array. Returned as void*, so it must be typecasted in use! */
void *Array_getDataPointer(Array a);

void Array_set_undef_value(Array a, double value); // readers may set the undef value
double Array_get_undef_value(Array a); 
/** check if array (after reading data from file!) element is the undefined value (also set by the reader) */
bool Array_is_data_undef(Array a, int idx);
/** check if whole array is undef.  Only after array is filled and calcminmax was executed!  */ 
bool Array_is_data_undef_all(Array a);
/** Generate an array with values from..to with stepping step */
Array Array_indgen(int from, int to, int step); 
/** Create transpose of a 2D array (allocates memory) */
Array Array_transpose(Array a);

void Array_zero(Array a);                         // fill with zeros

// debug
void Array_print(Array a); //print to screen

#endif
