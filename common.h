#ifndef __COMMON_H__
#define __COMMON_H__

#define myfree(p) if (p) { free(p); p=NULL; }

#ifndef HAVE_STRNDUP
#  define strndup(str,len) strdup(str)
#endif

#ifndef min
#  define min(x1,x2) ((x1) > (x2) ? (x2):(x1)) 
#endif
#ifndef max
#  define max(x1,x2) ((x1) > (x2) ? (x1):(x2))
#endif

#ifndef __cplusplus
    typedef char bool;  // good until sizeof(bool) in C++ = sizeof(char) in C
#endif
#define false 0
#define true  1

#define CUT_TO_BYTE(x) (x < 0 ? 0 : (x > 255 ? 255 : x)) 

typedef enum { undef_file, netcdf_file, hdf5_file, adiosbp_file, ascii_file } FileType;

#ifdef __MAIN__
char const * const fileTypeStrings[] = { "UNDEFINED", "NETCDF", "HDF5", "ADIOS", "ASCII" };
int verbose;
bool print_provenance;
#else
extern char const * const fileTypeStrings[];
extern int verbose;
extern bool print_provenance;
#endif


#define MAX_DIMS 5

#endif
