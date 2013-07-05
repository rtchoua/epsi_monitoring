#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#ifndef __USE_ISOC99
#  define __USE_ISOC99
#endif
#include <math.h>         // ISO_C99 defines INFINITY
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "common.h"

// global vars for this module
static FILE* fminmax = NULL;    // minmax file descriptor

#define GTABLE_MAXSIZE 1000

typedef struct {
    char *name;
    double gmin, gmax;
} gelem;

static gelem gtable[GTABLE_MAXSIZE];
static int elemnum;

// add new element or if already exists calculate new global min/max
static int gtable_add(char * name, double lmin, double lmax, double undef) {
    int i;
    for (i = 0; i < elemnum; i++) {
	if (!strcmp(name, gtable[i].name)) break;	
    }
    if (i == GTABLE_MAXSIZE) {
	fprintf(stderr, "Error: minmax: table size %d overflow. Recompile with larger table.\n", i);
	exit(33);
    }
    if (i == elemnum) {
	// add new item
	gtable[i].name = strdup(name);
	gtable[i].gmin = lmin;
	gtable[i].gmax = lmax;
	elemnum++;
    } else {
	// recalc min/max
	if (gtable[i].gmin == undef || lmin < gtable[i].gmin) gtable[i].gmin = lmin;
	if (gtable[i].gmax == undef || lmax > gtable[i].gmax) gtable[i].gmax = lmax;
    }
    return i;
}

// get element 
static int gtable_get(char * name) {
    int i;
    for (i = 0; i < elemnum; i++) {
	if (!strcmp(name, gtable[i].name)) break;	
    }
    if (i == elemnum) {
	// no such item
	i = -1;
    }
    return i;
}

static void init_gtable(void)
{
  memset(gtable, 0, GTABLE_MAXSIZE*sizeof(gelem));
  elemnum=0;
}




/**********************************************************/
/*                                                        */ 
/*                MIN / MAX FUNCTIONS                     */ 
/*                                                        */ 
/**********************************************************/
bool minmax_open( char *path, bool readmode) {  // mode: 1=read, 0=append 
    char fn[1024];
    if ( path != NULL) 
	strncpy(fn, path, 1024);
    else
	snprintf(fn, 1024, "minmax1D.dat");

    if ( readmode ) {
	if ((fminmax = fopen(fn,"r")) == NULL) {
	    fprintf(stderr, "Error at opening for reading minmax file %s: %s\n",
		    fn, strerror(errno));
	    return false;
	}
    } else {
	if ((fminmax = fopen(fn,"a+")) == NULL) {
	    fprintf(stderr, "Error at opening for writing minmax file %s: %s\n",
		    fn, strerror(errno));
	    return false;
	}
    }

    init_gtable();
    return true;
}

static void read_minmax(double undef) {
    int n, ts;
    double lmin, lmax;
    char name[1024];
    // read in the values into gtable
    if (fminmax) {
        fseek(fminmax, 0, SEEK_SET);
	while ( (n = fscanf( fminmax, "%[^;];%d;%lg;%lg\n", name, &ts, &lmin, &lmax)) != EOF) {
	    gtable_add(name, lmin, lmax, undef);
            //printf("Read local min/max for %s %d: %.12g   %.12g\n", name, ts, lmin, lmax);
	}
    }
}

void minmax_close() { 
    if (fminmax)
	fclose(fminmax); 
}

void minmax_print(char* varname, int time, double min, double max) {
    if (fminmax)
	fprintf(fminmax, "%s;%4.4d;%.12g;%.12g\n", varname, time, min, max);
}

void minmax_get(char* varname, double undef_value, double *gmin, double *gmax) {
    // get the global min/max value for a variable
    int n;
    char line[1024];
    double undef;
    //double UNDEF=9969209968386869046778552952102584320.000000;
    //double UNDEF=NC_FILL_DOUBLE;

    // get the undef value formatted as in the file
    sprintf(line, "%.12g", undef_value); 
    sscanf(line, "%lg", &undef);
    //printf("undef_value=%.12g,  str=%s, undef=%.12g\n", undef_value, line, undef);
    
    if (elemnum == 0) {
	read_minmax(undef);  // at first call, read in the minmax values from file
    }
    
    n = gtable_get(varname); 

    if (n > -1) {
	*gmin = gtable[n].gmin;
	*gmax = gtable[n].gmax;
	// convert back to full undef value if that was the case
	if (*gmin == undef) *gmin = undef_value;
	if (*gmax == undef) *gmax = undef_value;
    } else {
	fprintf(stderr, "Error: minmax_get: variable %s not found in minmax data\n", varname);
        *gmin = undef_value;
	*gmax = undef_value;
    }
    //if (verbose) printf("Global min/max for %s: %.13g   %.12g\n", varname, *gmin, *gmax);
}

