#ifndef __READADIOSBP_H__
#define __READADIOSBP_H__

#include "common.h"
#include "array.h"
#include "reader.h"



/** Open file and read info on all variables. 
  * Return DataFile. On error, DataFile.fd == -1
  */
DataFile readadiosbp_open(char *path, char *nameattr); 

/** Close data file and free all allocated resources (names and varinfos) */
void  readadiosbp_close(DataFile df);   


/** read the data of one variable with start[] and count[] hyperslabbing. 
  * The array will be initialized and allocated (so it should be freed manually after use).
  * On error: Array.type = undefArray;
  */
Array readadiosbp_read_var( DataFile df, VarInfo vi, int *start, int *count); 


#endif
