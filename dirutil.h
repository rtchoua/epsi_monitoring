#ifndef __DIRUTUL_H_
#define __DIRUTUL_H_

#include "common.h"

char* getbasename( char *path, char **dirname );

bool file_exists (char *path);       // true if stat(path) succeeds, i.e. it is an accessible item on file system
bool is_dir (char *path);            // true if path is a directory (and stat() succeeds on it)

/** mkdir -r
  * return: 0 on success, otherwise mkdir() syscall's return value
  */
int createdir_recursive( char* path);

#endif 
