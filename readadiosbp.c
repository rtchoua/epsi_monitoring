/* ADIOS BP data reader for Array
 *
 * Copyright Oak Ridge National Laboratory 2009
 * Author: Norbert Podhorszki, pnorbert@ornl.gov
**/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include "adios_read.h"
#include "adios_types.h"
#include "readadiosbp.h"
#include "dirutil.h"

/*
static void print_namelist (char **namelist, int length) {
    int i;
    for (i=0; i<length; i++)
        printf("\t%d: \t%s\n", i, namelist[i]);
    return;
}
*/

static void bpexit(int code, ADIOS_FILE *fp, ADIOS_GROUP *gp) {
    if (gp)
        adios_gclose (gp);
    if (fp)
        adios_fclose (fp);
    exit(code);
}


DataFile readadiosbp_open(char *path, char *nameattr) {
    ADIOS_FILE    * fp;
    ADIOS_GROUP   * gp;
    ADIOS_VARINFO * avi;
    int     gr, vr, i, j;             // loop vars
    int     mpi_comm_dummy;
    DataFile df = DataFile_new();
    
    if (verbose) printf("\nADIOS BP open: read header info from %s\n", path);

    // open the BP file
    // get number of groups, variables, timesteps, and attributes 
    fp = adios_fopen (path, mpi_comm_dummy); // only rank=0 reads in anything
    if (fp == NULL) {
	fprintf(stderr, "readadiosbp: error opening bp file %s:\n%s", path, adios_errmsg());
	bpexit(7, NULL, NULL);
    }

    if (verbose) {
         printf ("# of groups: %d\n", fp->groups_count); 
         printf ("# of variables: %d\n", fp->vars_count);
         printf ("# of attributes: %d\n", fp->attrs_count);
         printf ("time steps: %d steps from %d\n", fp->ntimesteps, fp->tidx_start);
    }

    // create VarInfo array for all variables
    df.fh         = (int64_t)fp;
    df.nvars      = fp->vars_count;
    df.ngroups    = fp->groups_count;
    df.time_start = fp->tidx_start-1;
    df.time_stop  = fp->tidx_start+fp->ntimesteps-1;
    df.varinfo    = VarInfo_newarray(fp->vars_count);
    if (df.varinfo == NULL) {
        fprintf(stderr, "error allocating varinfo array, nvars=%d\n", fp->vars_count);
        bpexit(7, fp, 0);
    }

    // each group has to be handled separately
    i = 0;  // var index for df.varinfo
    for (gr=0; gr<fp->groups_count; gr++) {
        if (verbose>1) printf("  group %s:\n", fp->group_namelist[gr]);
        gp = adios_gopen_byid (fp, gr);
        if (gp == NULL) {
	    fprintf(stderr, "readadiosbp: error opening group %s:\n%s", fp->group_namelist[gr], adios_errmsg());
	    bpexit(8, fp, NULL);
        }

        for (vr=0; vr<gp->vars_count; vr++) {
            avi = adios_inq_var_byid (gp, vr);
            if (avi == NULL) {
	        fprintf(stderr, "readadiosbp: error inquiring variable %s:\n%s", gp->var_namelist[vr], adios_errmsg());
	        bpexit(9, fp, gp);
            }
            //debug
            if (verbose>2) {
                printf("    v:   %s\tdims=%d (", gp->var_namelist[vr], avi->ndim);
                for (j=0; j<avi->ndim; j++) printf("%llu ", avi->dims[j]);
                printf(")\n");
            }
            // error check dimension size
            if (avi->ndim > MAX_DIMS) {
                fprintf(stderr, "error: variable %s/%s has more dimensions (%d) than supported (%d).\n", 
                        fp->group_namelist[gr], gp->var_namelist[vr], avi->ndim, MAX_DIMS);
                fprintf(stderr, "Recompile plotter with a larger reader.h:MAX_DIMS\n");
            }
            df.varinfo[i].name = strdup(gp->var_namelist[vr]);
            df.varinfo[i].dispname  = df.varinfo[i].name;  // default is full var path
            df.varinfo[i].groupname = strdup(fp->group_namelist[gr]);
            df.varinfo[i].varid = vr; 
            df.varinfo[i].ndims = min(avi->ndim, MAX_DIMS);
            df.varinfo[i].timedim = avi->timedim;

            // copy dimension sizes into varinfo.dimsize array
            if (verbose>1) printf("    var: %s\tdims=%d (", gp->var_namelist[vr], avi->ndim);
            for (j=0; j<df.varinfo[i].ndims; j++) { 
                df.varinfo[i].dimsize[j] = avi->dims[j];
                if (verbose>1) printf("%s%d ", (j==df.varinfo[i].timedim ? "T-": ""), df.varinfo[i].dimsize[j]);
            }
            if (verbose>1) printf(")\n");

            // determine the type of the variable
            switch(avi->type) {
                case adios_unsigned_byte:  
                case adios_byte:
                case adios_string:
                    df.varinfo[i].type = int8Array;
                    break;
                        
                case adios_unsigned_short:  
                case adios_short:
                    df.varinfo[i].type = int16Array;
                    break;
                        
                case adios_unsigned_integer:
                case adios_integer:    
                        df.varinfo[i].type = int32Array;
                    break;

                case adios_unsigned_long:
                case adios_long:        
                        df.varinfo[i].type = int64Array;
                    break;

                case adios_real:
                    df.varinfo[i].type = floatArray;
                    break;

                case adios_double:
                case adios_long_double:
                    df.varinfo[i].type = doubleArray;
                    break;

                case adios_complex:  
                case adios_double_complex:
                default:
                    fprintf(stderr, "Adios type %d (%s) not supported in plotter. group=%s var=%s\n", 
                        avi->type, adios_type_to_string(avi->type), fp->group_namelist[gr], gp->var_namelist[vr]);
                    df.varinfo[i].type = undefArray;
                    break;
            }

            i++; // next varinfo
            adios_free_varinfo(avi);
        }
        adios_gclose (gp);
    }

    return df;
}                

/** Close data file and free all allocated resources (names and varinfos) */
void  readadiosbp_close(DataFile df) {
    adios_fclose ( (ADIOS_FILE *)df.fh);
}   


/** read the data of one variable of one timestep.
  * The array will be initialized and allocated (so it should be freed manually after use).
  * On error: Array.type = undefArray;
*/
Array readadiosbp_read_var( DataFile df, VarInfo vi, int *start, int *count) {
    Array x = Array_new(undefArray);
    int i,j;
    ADIOS_GROUP *gp;
    ADIOS_VARINFO *avi;
    // to check data type
    ArrayType atype;
    // to read data
    uint64_t start_t[MAX_DIMS], count_t[MAX_DIMS];
    int nelems;   // number of elements for one read
    int elemsize; // size in bytes of one number
    int st, ct;
    void *data;
    char buf[64];
    int  ndims, dims[MAX_DIMS]; // resulting array dimensions
    uint64_t bytes_read;

    // open the group 
    if (verbose>2) printf("  read group=%s var=%s\n", vi.groupname, vi.name);
    gp = adios_gopen ( (ADIOS_FILE *)df.fh, vi.groupname);
    
    // get type of var again from the bp file directly 
    // and set the array accordingly
    avi = adios_inq_var_byid (gp, vi.varid);
    switch(avi->type) {
        case adios_unsigned_byte:
        case adios_byte:
        case adios_string:
             atype = int8Array;
             break;

        case adios_unsigned_short:  
        case adios_short:
            atype = int16Array;
            break;
                
        case adios_unsigned_integer:
        case adios_integer:    
            atype = int32Array;
            break;

        case adios_unsigned_long:
        case adios_long:        
            atype = int64Array;
            break;

        case adios_real:
            atype = floatArray;
            break;

        case adios_double:
        case adios_long_double:
            atype = doubleArray;
            break;

        case adios_complex:  
        case adios_double_complex:
        default:
            fprintf(stderr, "Adios type %d (%s) not supported in plotter. var=%s\n", 
                    avi->type, adios_type_to_string(avi->type), vi.name);
	    Array_seterrno(x,-3);
            adios_gclose(gp);
	    return x;
    }
    Array_setType(x, atype);
    elemsize = arrayTypeElementSizes[ atype ]; // 1-2-4-8 bytes

    // debug 
    if(verbose>2) printf("-- debug: type of array x should be %s\n", arrayTypeStrings[atype]);
    if(verbose>2) printf("-- debug: type of array x is %s\n", arrayTypeStrings[Array_getType(x)]);

    // create the counter arrays with the appropriate lengths
    // transfer start and count arrays to format dependent arrays
    nelems  = 1;
    ndims = 0;
    for (i=0; i<vi.ndims; i++) {
        if (start[i] < 0)  // negative index means last-|index|
            st = vi.dimsize[i]+start[i];
        else
            st = start[i];
        if (count[i] < 0)  // negative index means last-|index|+1-start
            ct = vi.dimsize[i]+count[i]+1-st;
        else
            ct = count[i];

        if (verbose>2) 
            printf("    i=%d, st=%d ct=%d\n", i, st, ct);

        start_t[i] = st;
        count_t[i] = ct;
        nelems *= ct;
        if (verbose>1) 
            printf("    s[%d]=%llu, c[%d]=%llu, n=%d\n", i, start_t[i], i, count_t[i], nelems);
        if (verbose>2) 
            printf("    i=%d, start_t[%d]=%lld count_t[%d]=%lld\n", i, i, start_t[i], i, count_t[i]);
        // build output array dimensions
        if (count_t[i] > 1) { 
            dims[ndims] = count_t[i];
            ndims++;
        }
    }

    if (verbose>1) {
        printf("    resulting array has %d dimensions: %d", ndims, dims[0]);
        for (i=1; i<ndims; i++) printf( "x%d", dims[i] );
        printf(" total size=%d\n", nelems);
    }
    
    // allocate array and get hold to it raw
    if (!Array_allocate_dims(x, ndims, dims)) {
	fprintf(stderr, "Error: readadiosbp_read_var: could not allocate array of %d elements when reading variable %s\n", 
                (int)nelems, vi.name);
	Array_seterrno(x,-5);
        adios_gclose(gp);
	return x;
    }
    //Array_zero(x);
    data = Array_getDataPointer(x);

    if (verbose>2) printf("adios_get_var name=%s data=%ld\n", vi.name, (long) data);
    bytes_read = adios_read_var_byid (gp, vi.varid, start_t, count_t, data);
    if (bytes_read < 0) {
        fprintf(stderr,"Error reading variable %s: %s\n", vi.name, adios_errmsg());
    } else {
        if (verbose>1) {
            printf("  Has read in %lld bytes\n", bytes_read);
        } 
        if (verbose>2) {
            printf("  data = [");
            for(j=0; j<nelems; j++) {
                Array_strvalue(x, j, buf);
                printf("%s ", buf);
            }
            printf("]\n");
        }
    }

    // FIXME: need error checking from adios !
    /*
    if (status != NC_NOERR) {
	fprintf(stderr, "error reading data for variable %d %s: %s\n",
	       vi.varid, vi.name, nc_strerror(status));
	Array_seterrno(x,-6);
	Array_free(x);
	return x;
    } */

    //FIXME: set undef value of ADIOS arrays
    Array_set_undef_value(x, 0.0);
    
    adios_free_varinfo(avi);
    adios_gclose(gp);
    return x;
}

