/*
 * Copyright Oak Ridge National Laboratory 2009
 * Authors: Norbert Podhorszki, pnorbert@ornl.gov
**/
#define __ARRAY_MAIN_

#include "common.h"
#include "array.h"
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> // memset
#ifndef __USE_ISOC99
#  define __USE_ISOC99
#endif
#include <math.h>

#define EXIT_ON_NULL(a, funcname)  if (!a) { \
        fprintf(stderr, "%s: called with NULL\n", funcname); \
        exit(110); \
    } 

char const * const arrayTypeStrings[] = { "undef", "int8", "int16", "int32", "int64", "float", "double" };
size_t const arrayTypeElementSizes[] = { 0, sizeof(char), sizeof(int16_t), sizeof(int32_t), sizeof(int64_t), sizeof(float), sizeof(double) };

struct ArrayStruct {
    ArrayType type;
    long      length;            // = dimsize[0] * ... * dimsize[ndims]
    int       ndims;             // number of dimensions
    int       dimsize[MAX_DIMS]; // size of each dimension
    double    min;
    double    max;
    double    undef_value;       // The undef value of the data set (should set by readers)
    int       err;               // can be set to non-zero by programs to indicate error
    union {
        struct {
            char *data;
        } int8;
        struct {
            int16_t *data;
        } int16;
        struct {
            int32_t *data;
        } int32;
        struct {
            int64_t *data;
        } int64;
        struct {
            float   *data;
        } realsp;
        struct {
            double  *data;
        } realdp;
    };
};

Array Array_transpose(Array a) {
    Array aT;
    bool b;
    int i,j, src, target;
    size_t esize;
    void *aTdata, *adata;

    EXIT_ON_NULL(a, "Array_transpose")

    aT = Array_new( a->type );

    if (a->ndims != 2) {
        fprintf(stderr, "Only 2D arrays can be transposed\n");
        return aT;
    }

    b = Array_allocate(aT, a->length);
    if (!b) {
        fprintf(stderr, "Could not allocate space for transposed array");
        return aT;
    }

    // flip dimensions
    aT->ndims = a->ndims;
    for (i=0; i < a->ndims; i++) {
        aT->dimsize[i] = a->dimsize[a->ndims-i-1];
    }

    // copy attributes
    aT->length = a->length;
    aT->min = a->min;
    aT->max = a->max;
    aT->undef_value = a->undef_value;
    
    // transpose data
    aTdata = Array_getDataPointer(aT);
    adata  = Array_getDataPointer(a);
    esize = arrayTypeElementSizes[a->type];
    target = 0;
    //printf("%%%%%%%%%%%%%%%%%% TRANSPOSE %%%%%%%%%%%%%%%%%%\n");
    //printf("esize = %d\n", esize);
    for (i=0; i < aT->dimsize[0]; i++) {
        for (j=0; j < aT->dimsize[1]; j++) {
            src = j*aT->dimsize[0] + i;
            // aT[i,j] = a[j,i] = a( j*M + i ), a=NxM, aT=MxN
            //printf("%d: src = %d\n", target, src);
            memcpy( aTdata + target*esize, adata + src*esize, esize);
            target++;
        }
    }
    return aT;
}

Array Array_new( ArrayType type) {
    Array a = (Array) malloc(sizeof(struct ArrayStruct));
    if (a) {
        a->type = type;
        a->length = 0;
        a->ndims  = 0;
        a->min = INFINITY;
        a->max = -INFINITY;
        a->undef_value = INFINITY;  // this is not valid, should be set by the data readers
        a->err = 0;
   }
   return a;
}

void Array_setType(Array a, ArrayType type) { 
    EXIT_ON_NULL(a, "Array_setType")
    a->type = type; 
}

ArrayType Array_getType(Array a) { 
    EXIT_ON_NULL(a, "Array_getType")
    return a->type; 
}

int Array_errno(Array a) {
    EXIT_ON_NULL(a, "Array_errno")
    return a->err;
}

void Array_seterrno(Array a, int err) {
    EXIT_ON_NULL(a, "Array_seterrno")
    a->err = err;
}

bool Array_allocate_dims(Array a, int ndims, int *dims) {
    EXIT_ON_NULL(a, "Array_allocate_dims")
    bool retval = false;
    int length = 1;
    int i;
    for (i=0; i<ndims; i++) length = length * dims[i];
    retval = Array_allocate(a, length);
    if (retval)
        retval = Array_setDimensions(a, ndims, dims);
    /*if (verbose>2) {
        printf("    Allocated array of size %d, dimensions: %d", length, a->dimsize[0]);
        for(i=1; i<ndims; i++) printf("x%d", a->dimsize[i]);
        printf("\n");
    }*/
    return retval;
}


bool Array_setDimensions(Array a, int ndims, int *dims) {
    EXIT_ON_NULL(a, "Array_setDimensions")
    bool retval = false;
    int length = 1;
    int i;
    for (i=0; i<ndims; i++) length = length * dims[i];
    // dims[0]*...*dims[ndims] must = Array_length(a)
    if (a->length == length) {
        a->ndims = ndims;
        for (i=0; i<ndims; i++) a->dimsize[i] = dims[i];
        retval = true;
    } else {
        fprintf(stderr, "Array_setDimensions: Length of array does not equal the productum of given dimensions\n");
    }
    return retval;
}

int Array_getDimensions(Array a) {
    // return ndims, and fill dims[] with the values. no alloc inside!
    EXIT_ON_NULL(a, "Array_getDimensions")
    return a->ndims;
}

int Array_getAllDimensions(Array a, int *dims) {
    // return ndims, and fill dims[] with the values. no alloc inside!
    EXIT_ON_NULL(a, "Array_getAllDimensions")
    int i;
    if (dims != NULL) 
        for (i=0; i<a->ndims; i++) 
            dims[i] = a->dimsize[i];
    return a->ndims;
}

int Array_getDimensionSize(Array a, int dimidx) {
    // return dimsize[dimidx] or -1 if idx is out of scope
    EXIT_ON_NULL(a, "Array_getDimensionSize")
    if (dimidx >= 0 && dimidx <= a->ndims) 
        return a->dimsize[dimidx];
    else {
        fprintf(stderr, "Array_getDimension: index %d is out of array dimensions 0..%d\n", dimidx, a->ndims);
        return -1;
    }
}

bool Array_allocate(Array a, int n) {
    EXIT_ON_NULL(a, "Array_allocate")
    bool retval = false;
    switch (a->type) {
        case int8Array:
            a->int8.data = malloc( n*sizeof(char) );
            retval = (a->int8.data != NULL);
            break;
        case int16Array:
            a->int16.data = malloc( n*sizeof(int16_t) );
            retval = (a->int16.data != NULL);
            break;
        case int32Array:
            a->int32.data = malloc( n*sizeof(int32_t) );
            retval = (a->int32.data != NULL);
            break;
        case int64Array:
            a->int64.data = malloc( n*sizeof(int64_t) );
            retval = (a->int64.data != NULL);
            break;
        case floatArray:
            a->realsp.data = malloc( n*sizeof(float) );
            retval = (a->realsp.data != NULL);
            break;
        case doubleArray:
            a->realdp.data = malloc( n*sizeof(double) );
            retval = (a->realdp.data != NULL);
            break;
        default:
            fprintf(stderr, "Array: Error: unknown array type used at Array_allocate\n");
            exit(100);
    }
    if (retval) {
        a->length = n;
        a->ndims  = 1;
        a->dimsize[0] = n;
        a->err = 0;
    }
    if (!retval) 
        a->err = -1;
    return retval;
}

void Array_free(Array a) {
    if (a) {
        switch (a->type) {
           case int8Array:
                myfree(a->int8.data);
                break;
           case int16Array:
                myfree(a->int16.data);
                break;
            case int32Array:
                myfree(a->int32.data);
                break;
            case int64Array:
                myfree(a->int64.data);
                break;
            case floatArray:
                myfree(a->realsp.data);
                break;
            case doubleArray:
                myfree(a->realdp.data);
                break;
            default:
                fprintf(stderr, "Array: Error: unknown array type used at Array_free\n");
                exit(101);
        }
        a->length = 0;
        a->ndims = 0;
        a->err = 0;
    }
    myfree(a);
}

long Array_length(Array a) { 
    EXIT_ON_NULL(a, "Array_length")
    return a->length; 
}


#define MINMAX(a, arr)   for (i = 0; i < a->length; i++) {    \
                             if (a->arr[i] < a->min) a->min = a->arr[i]; \
                             if (a->arr[i] > a->max) a->max = a->arr[i]; \
                         }

void Array_calcminmax(Array a) {
    long i;
    EXIT_ON_NULL(a, "Array_calcminmax")
    switch (a->type) {
        case int8Array:
            MINMAX(a, int8.data)
            break;
        case int16Array:
            MINMAX(a, int16.data)
            break;
        case int32Array:
            MINMAX(a, int32.data)
            break;
        case int64Array:
            MINMAX(a, int64.data)
            break;
        case floatArray:
            MINMAX(a, realsp.data)
            break;
        case doubleArray:
            MINMAX(a, realdp.data)
            break;
        default:
            fprintf(stderr, "Array: Error: unknown array type used at Array_free\n");
            exit(102);
    }
    //printf("DEBUG: Array_calcminmax: min = %le, max = %le undef = %le\n", 
    //            (double) a->min, (double) a->max, a->undef_value);
}


double Array_min(Array a) {
    EXIT_ON_NULL(a, "Array_min")
    if (a->min == INFINITY && a->max == -INFINITY) {
        //fprintf(stderr, "Array: Array_min should be called after Array_calcminmax only\n");
        Array_calcminmax(a);
    }
    return a->min;

}

double Array_max(Array a) {
    EXIT_ON_NULL(a, "Array_max")
    if (a->min == INFINITY && a->max == -INFINITY) {
        fprintf(stderr, "Array: Array_max should be called after Array_calcminmax only\n");
        Array_calcminmax(a);
    }
    return a->max;
}

double Array_doublevalue(Array a, long idx) {
    EXIT_ON_NULL(a, "Array_doublevalue")
    switch (a->type) {
        case int8Array:
            return (double) a->int8.data[idx];
            break;
        case int16Array:
            return (double) a->int16.data[idx];
            break;
        case int32Array:
            return (double) a->int32.data[idx];
            break;
        case int64Array:
            return (double) a->int64.data[idx];
            break;
        case floatArray:
            return (double) a->realsp.data[idx];
            break;
        case doubleArray:
            return a->realdp.data[idx];
            break;
        default:
            fprintf(stderr, "Array: Error: unknown array type used at Array_free\n");
            exit(102);
    }

}

void Array_strvalue(Array a, long idx, char* str) {
    EXIT_ON_NULL(a, "Array_strvalue")
    switch (a->type) {
        case int8Array:
            sprintf(str, "%d", a->int8.data[idx]);
            break;
        case int16Array:
            sprintf(str, "%d", a->int16.data[idx]);
            break;
        case int32Array:
            sprintf(str, "%d", a->int32.data[idx]);
            break;
        case int64Array:
            sprintf(str, "%lld", (long long int)a->int64.data[idx]);
            break;
        case floatArray:
            sprintf(str, "%.15lf", a->realsp.data[idx]);
            break;
        case doubleArray:
            sprintf(str, "%.15le", a->realdp.data[idx]);
            break;
        default:
            fprintf(stderr, "Array: Error: unknown array type used at Array_free\n");
            exit(102);
    }
}

void Array_strvalue_padded(Array a, long idx, char* str, int width) {
    EXIT_ON_NULL(a, "Array_strvalue")
    switch (a->type) {
        case int8Array:
            sprintf(str, "%0*d", width, a->int8.data[idx]);
            break;
        case int16Array:
            sprintf(str, "%0*d", width, a->int16.data[idx]);
            break;
        case int32Array:
            sprintf(str, "%0*d", width, a->int32.data[idx]);
            break;
        case int64Array:
            sprintf(str, "%0*lld", width, (long long int)a->int64.data[idx]);
            break;
        case floatArray:
            sprintf(str, "%0*lf", width, a->realsp.data[idx]);
            break;
        case doubleArray:
            sprintf(str, "%0*le", width, a->realdp.data[idx]);
            break;
        default:
            fprintf(stderr, "Array: Error: unknown array type used at Array_free\n");
            exit(102);
    }
}

void *Array_getDataPointer(Array a) {
    void *data;
    EXIT_ON_NULL(a, "Array_getDataPointer")
    switch (a->type) {
        case int8Array:
            data = (void *) a->int8.data;
            break;
        case int16Array:
            data = (void *) a->int16.data;
            break;
        case int32Array:
            data = (void *) a->int32.data;
            break;
        case int64Array:
            data = (void *) a->int64.data;
            break;
        case floatArray:
            data = (void *) a->realsp.data;
            break;
        case doubleArray:
            data = (void *) a->realdp.data;
            break;
        default:
            fprintf(stderr, "Array: Error: unknown array type used at Array_free\n");
            exit(102);
    }
    return data;
}
 
void Array_set_undef_value(Array a, double value) {
    EXIT_ON_NULL(a, "Array_set_undef_value")
    a->undef_value = value;
}

double Array_get_undef_value(Array a) {
    EXIT_ON_NULL(a, "Array_get_undef_value")
    return a->undef_value;
}


/** check if array (after reading data from file!) element is the undefined value (also set by the reader)
  */
bool Array_is_data_undef(Array a, int idx) {
    bool retval = false;
    EXIT_ON_NULL(a, "Array_is_data_undef")
    switch (a->type) {
        case int8Array:
            retval = ((double)a->int8.data[idx] == a->undef_value);
            break;
        case int16Array:
            retval = ((double)a->int16.data[idx] == a->undef_value);
            break;
        case int32Array:
            retval = ((double)a->int32.data[idx] == a->undef_value);
            break;
        case int64Array:
            retval = ((double)a->int64.data[idx] == a->undef_value);
            break;
        case floatArray:
            retval = ((double)a->realsp.data[idx] == a->undef_value);
            break;
        case doubleArray:
            retval = ((double)a->realdp.data[idx] == a->undef_value);
            break;
        default:
            fprintf(stderr, "Array: Error: unknown array type used at Array_free\n");
            exit(102);
    }
    return retval;
}

/** check if whole array is undef.
    Only after array is filled and calcminmax was executed!
  */
bool Array_is_data_undef_all(Array a) {
    EXIT_ON_NULL(a, "Array_is_data_undef_all")
    //return ((double)a->min == a->undef_value && (double)a->max == a->undef_value);
    return (Array_min(a) == a->undef_value && Array_max(a) == a->undef_value);
    // Note: Array_min(a) will force calculating min if not yet calculated !
}

Array Array_indgen(int from, int to, int step) {
    int i=0, v;
    int dims[1];
    Array a = Array_new(int32Array);
    //printf("indgen: [%d..%d] with step %d. Type = %d\n", from, to, step, Array_getType(a));
    if (step < 0) {
        fprintf(stderr, "Array_indgen: negative stepping is not supported\n");
        exit(21);
    }
    if (from > to) {
        fprintf(stderr, "Array_indgen: from (%d) > to (%d)\n", from, to);
        exit(20);
    }
    if (a) {
        if (Array_allocate(a, (to - from) / step + 1)) {
            for (v=from; v<=to; v+=step) {
                a->int32.data[i] = v;
                i++;
            }
            a->min = (double) from;
            a->max = (double) a->int32.data[(to - from) / step];
            //printf("indgen: array created, len=%d (%d)\n", Array_length(a), a->length);
        } else {
            fprintf(stderr, "Array_indgen: cannot allocate array for %d int32 elements\n", (to - from) / step + 1);
        }
    }
    dims[0] =  (to - from) / step + 1;
    Array_setDimensions(a, 1, dims);
    //printf("indgen: array ptr = %x\n", a);
    return a;
}

void Array_zero(Array a) {
    EXIT_ON_NULL(a, "Array_zero")
    switch (a->type) {
        case int8Array:
            memset(a->int8.data, 0, (size_t)Array_length(a) * sizeof(char));
            break;
        case int16Array:
            memset(a->int16.data, 0, (size_t)Array_length(a) * sizeof(int16_t));
            break;
        case int32Array:
            memset(a->int32.data, 0, (size_t)Array_length(a) * sizeof(int32_t));
            break;
        case int64Array:
            memset(a->int64.data, 0, (size_t)Array_length(a) * sizeof(int64_t));
            break;
        case floatArray:
            memset(a->realsp.data, 0, (size_t)Array_length(a) * sizeof(float));
            break;
        case doubleArray:
            memset(a->realdp.data, 0, (size_t)Array_length(a) * sizeof(double));
            break;
        default:
            fprintf(stderr, "Array: Error: unknown array type used at Array_free\n");
            exit(102);
    }

}

void Array_print(Array a) {
    int i;
    EXIT_ON_NULL(a, "Array_print")
    printf("===========================\nArray dump: lenght=%ld\n", a->length);
    for (i = 0; i < a->length; i++) {
        printf("%5d: %lf\n", i, Array_doublevalue(a,i)); 
    }
    printf("===========================\n");
}
