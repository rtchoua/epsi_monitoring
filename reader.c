/**
 * Copyright Oak Ridge National Laboratory 2009
 * Author: Norbert Podhorszki, pnorbert@ornl.gov
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "reader.h"  // reads common.h and array.h

#ifdef NETCDF
#    include "readcdf.h"    // netcdf reader
#endif
#ifdef ADIOS
#    include "readadiosbp.h"  // adios BP reader
#endif
#ifdef HDF5
#    include "readhdf5.h"   // hdf5   reader
#endif

/** Get from string the file type, i.e. netcdf or hdf5 or adios
  */
FileType reader_getFileType(char *str) {
    if (!strcasecmp(str,"netcdf") ||
        !strcasecmp(str,"cdf")    ||
        !strcasecmp(str,"nc") 
       )
        return netcdf_file;
    else if (!strcasecmp(str,"hdf5") ||
             !strcasecmp(str,"h5")
            )
        return hdf5_file;
    else if (!strcasecmp(str,"adiosbp") ||
             !strcasecmp(str,"bp")      ||
             !strcasecmp(str,"adios")
            )
        return adiosbp_file;
    else
        return undef_file;
}

// try to determine the file's type 
static FileType  reader_guessFileType(char *path) {
    FileType type = undef_file;
    char *ext = strrchr(path, (int)'.');
    if (verbose>1) printf("guess file type from ext: %s\n", ext);
    if (ext != NULL) {
       ext++; // get rid of '.'
       if (ext != NULL)
           type = reader_getFileType(ext);
    }
    return type;
}


VarInfo VarInfo_new(void) {
    VarInfo vi;
    vi.name = NULL;
    vi.dispname = NULL;
    vi.varid = -1;
    vi.ndims = 0;
    vi.timedim = -1;
    vi.groupname = NULL;
    return vi;
}

VarInfo* VarInfo_newarray(int n) {
    VarInfo *va;
    int i;
    if (n < 1) {
        fprintf(stderr,"Error: VarInfo_newarray() is called with an argument < 1: %d\n", n);
        return NULL;
    }
    va = malloc( n * sizeof(VarInfo) );
    if (va == NULL) {
        fprintf(stderr,"Error: Could not allocate %d VarInfo structs\n", n);
        return va;
    }
    for (i=0; i<n; i++) {
        va[i] = VarInfo_new();
    }
    return va;
}
    
DataFile reader_open(char *path, FileType type, char *nameattr) {
    DataFile df;
    if (type == undef_file) {
        printf("Make sure i'm not guessing\n");
	// Try to find out what is its type
	type = reader_guessFileType(path);
	if (type == undef_file) {
	    fprintf(stderr, "reader_open: could not determine file type for %s\n", path);
	    exit(200);
	}
    }
    switch (type) {
	case(netcdf_file):
#ifdef NETCDF
	    df = readcdf_open(path, nameattr);
	    df.type = netcdf_file;
#else
	    fprintf(stderr, "reader_open: NetCDF file type support was not compiled into plotter.\n");
	    exit(200);
#endif
	    break;
	case(adiosbp_file):
#ifdef ADIOS
	    df = readadiosbp_open(path, nameattr);
	    df.type = adiosbp_file;
#else
	    fprintf(stderr, "reader_open: ADIOS BP file type support was not compiled into plotter.\n");
	    exit(200);
#endif
	    break;
	case(hdf5_file):
#ifdef HDF5
	    df = readhdf5_open(path, nameattr);
	    df.type = hdf5_file;
#else
	    fprintf(stderr, "reader_open: HDF5 file type support was not compiled into plotter.\n");
	    exit(200);
#endif
	    break;
	case(ascii_file):
	    //df = readascii_open(path);
	    df.type = ascii_file;
	    break;
	default:
	    fprintf(stderr, "reader_open: unknown file type is given for %s : %d\n", path, (int) type);
	    exit(200);
	    break;
    }
    return df;
}

void reader_close(DataFile df) {
    int i;

    if (df.fd == -1 && df.fh == -1)
        return;

    switch (df.type) {
#ifdef NETCDF
	case(netcdf_file):
	    readcdf_close(df);
	    break;
#endif
#ifdef ADIOS
	case(adiosbp_file):
	    readadiosbp_close(df);
	    break;
#endif
#ifdef HDF5
	case(hdf5_file):
	    readhdf5_close(df);
	    break;
#endif
	case(ascii_file):
	    //readascii_close(df);
	    break;
	default:
	    fprintf(stderr, "reader_close: unknown file type is given: %d\n", (int) df.type);
	    exit(200);
	    break;
    }

    // free up strings in each var info
    for (i=0; i < df.nvars; i++) {
        if (df.varinfo[i].dispname != df.varinfo[i].name)
            myfree(df.varinfo[i].dispname);
	myfree(df.varinfo[i].name);
        df.varinfo[i].dispname = NULL; // make sure even if it was not freed (=name)
	myfree(df.varinfo[i].groupname);
    }

    // free up varinfo 
    myfree(df.varinfo);

    df.fd = -1;
    df.fh = -1;
}

VarInfo reader_getVarInfo( DataFile df, char *varname ) {
    int i = 0;
    int varstartpos = 0; // from where to match varname (beginning or +1)
    int dfstartpos = 0;  // from where to match df.varinfo[i] (beginning or +1)
    /* Find the variable: full path is stored with a starting / 
       We need to match names given with or without the starting /
       startpos is 0 or 1 to indicate if the argument has starting / or not
    */
    if (varname[0] == '/') varstartpos = 1;
    VarInfo vi = VarInfo_new();
    for (i=0; i < df.nvars; i++) {
        if (df.varinfo[i].name[0] == '/') dfstartpos = 1;
        //printf("getVarInfo: varname=%s, df.varinfo[%d]=%s, vstart=%d, dfstart=%d",
        //        varname,  df.varinfo[i].name, varstartpos, dfstartpos);
	if (!strcmp(df.varinfo[i].name+dfstartpos, varname+varstartpos))
	    return df.varinfo[i];
    }
    // not found, return a dummy
    return vi; // vi.name == NULL
}


Array reader_read_var( DataFile df, VarInfo vi, int *start, int *count)  {
    Array a;
    switch (df.type) {
#ifdef NETCDF
	case(netcdf_file):
	    a = readcdf_read_var(df, vi, start, count);
	    break;
#endif
#ifdef ADIOS
	case(adiosbp_file):
	    a = readadiosbp_read_var(df, vi, start, count);
	    break;
#endif
#ifdef HDF5
	case(hdf5_file):
	    a = readhdf5_read_var(df, vi, start, count);
	    break;
#endif
	case(ascii_file):
	    //a = readascii_read_var(df, vi, start, count);
	    break;
	default:
	    fprintf(stderr, "reader_close: unknown file type is given: %d\n", (int) df.type);
	    exit(200);
	    break;
    }
    return a;
}

Array reader_read_var_byname( DataFile df, char *varname, int *start, int *count)  {
    Array a;

    int vidx;
    // look for var in the varinfo array
    for (vidx = 0; vidx < df.nvars; vidx++) {
        if (!strcmp(varname, df.varinfo[vidx].name))
            break;
    }

    if (vidx < df.nvars) {
        a = reader_read_var(df, df.varinfo[vidx], start, count);
    } else {
        // if not found, return with error
        fprintf(stderr, "Error: readhdf5_read_var: Variable %s not found in the file\n", varname);
        a = Array_new(undefArray);
        Array_seterrno(a,-1);
    }
    return a;
}

void print_DFinfo( DataFile df) {
    int i, j;
    printf("==============================\n"); 
    printf("DataFile: fd=%d", df.fd);
    printf(" type=%s", fileTypeStrings[df.type]);
    //printf(" timesteps=%d", df.timeTo-df.timeFrom+1);
    printf(" nvars=%d\n", df.nvars); 
    //printf("DataFile: fd=%d  type=%s  timesteps=%d  nvars=%d\n", 
    //	    df.fd, fileTypeStrings[df.type], df.timesteps, df.nvars); 
    for (i=0; i < df.nvars; i++) {
	printf("  var=%s  type=%d  varid=%d  ndims=%d (", 
		df.varinfo[i].name, df.varinfo[i].type, 
		df.varinfo[i].varid, df.varinfo[i].ndims);
	for (j=0; j < df.varinfo[i].ndims; j++) 
	    printf("%d ", df.varinfo[i].dimsize[j]);
	printf(")\n");
    }
    printf("==============================\n"); 
}

