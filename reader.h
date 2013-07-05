#ifndef __READER_H__
#define __READER_H__

#include <stdint.h>
#include "common.h"
#include "array.h"

#define MAX_DIMS 5

typedef struct {
    char *name;                  // name of var (with full path if available)
    char *dispname;              // name to use for display and minmax (=name or an attribute from the file)
    ArrayType type;              // type of var
    int varid;                   // integer id of var in datafile if applicable for the reader
    int ndims;                   // number of dimensions
    int dimsize[MAX_DIMS];       // size of each dimension
    int timedim;                 // ADIOS: 0:index of the time dimension for this variable, -1 if there is no time dim
    char *groupname;             // ADIOS: adios group of variable is not part of path, we need to store this info
} VarInfo;

typedef struct {
    // if you change order, or add a field then change the DataFile_new() macro too!
    FileType type;          //    = undef_file; 
    int      fd;            //    = -1;     // file descriptor
    int64_t  fh;            //    = -1;     // ADIOS only: file handler is 64 bit integer
    int      nvars;         //    = 0;      // number of variables (infos)
    int      ngroups;       //    = -1;     // ADIOS only: number of groups (to know if it is 1 or more)
    int      time_start;    //    = 0;      // ADIOS only: first time index in the file, 0 is not valid in ADIOS
                                            //    0 is the default because it is used in calculations for all file types
    int      time_stop;     //    = 0;      // ADIOS only: last time index in the file
    VarInfo  *varinfo;      //    = NULL;   // some info about each var
} DataFile;

/** to initialize a DataFile variable, instead of 
        DataFile df; 
    use
       	DataFile df = DataFile_new();
*/
#define DataFile_new() { undef_file, -1, -1, 0, -1, 0, 0, NULL} 

/** to initialize a VarInfo variable, instead of 
       VarInfo vi;
    use
       VarInfo vi = VarInfo_new();
  */
VarInfo VarInfo_new(void);

/** allocate an array of 'n' VarInfo structs initialized.
    uses malloc(), so use free() when not needed anymore
    returns NULL if allocation fails
  */
VarInfo* VarInfo_newarray(int n);

/** Open file and read info on all variables. 
  * If type is undef, we try to guess what it is.
  *
  * nameattr is an optional string which tells the readers that
  * an attribute of each variable object with nameattr should be 
  * looked up and used as variable display name instead of the
  * path of the variable. Use NULL if there is no such attribute.
  *
  * Return DataFile. On error, DataFile.fd == -1
  */
DataFile  reader_open(char *path, FileType type, char *nameattr); 

/** Close data file and free all allocated resources (names and varinfos) */
void      reader_close(DataFile df);   

/** Set the start[] and count[] arrays to read variables.
  * start:  "[t][-]n[,...]"; 't' in a dimension denotes the time dimension
  *                          only one t is allowed
  * count:  "[-]n[,...]"; -1 denotes 'until end'
  */
void reader_setCounters(DataFile df, char *start, char *count);

/** Get info about one variable.
  * If variable is not found: VarInfo.name == NULL
  */
VarInfo reader_getVarInfo( DataFile df, char *varname );  

/** read the data of one variable of one timestep.
  * The array will be initialized and allocated (so it should be freed manually after use).
  * On error: Array.type = undefArray;
  */
Array reader_read_var( DataFile df, VarInfo vi, int *start, int *count); 
Array reader_read_var_byname( DataFile df, char *varname, int *start, int *count); 

/********** PUBLIC but not necessary to use ******************/
FileType  reader_getFileType(char *str);   // get file type from the content of string
void print_DFinfo( DataFile df); // debug

#endif
