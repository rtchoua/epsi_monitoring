#ifndef _GNU_SOURCE
#   define _GNU_SOURCE
#endif

#define __MAIN__
#include "common.h"
#include "array.h"
#include "dirutil.h"
#include "graceplot.h"
#include "reader.h"
#include "minmax.h"
#include "timer.h"

#ifdef VTK
#   include "vtk-graph.h"
#endif

#include "adios_read.h"
#include "dirutil.h"
#include <dirent.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <errno.h>
#include <limits.h>   
#include <math.h>    
#include <libgen.h>   
#include <regex.h>   
#include <fnmatch.h>  
#include "adios_read.h"

// Global variables (for main.c only)
// Values from the arguments or defaults
char *outpath;              // output files' starting path (can be extended with subdirs, names, indexes)
char *inpath;               // input files containing data files on titan
char *shotname;             // shotname to be part of output dir on sith
char *xvar;                 // name of x variable
char *xvar1;                // name of x variable
char *tvar;                 // name of time variable
char *xfile;                // file containing x variable
char *tfile;                // file containing t variable
char *vfile;                // file containing variable (values) to plot
char *minmaxfile;           // file containing min/max values for variables (written by plotter)
char *start;                // dimension spec starting points
char *count;                // dimension spec counts
char *xstart;               // dimension spec starting points
char *xcount;               // dimension spec counts
char *nameattr;             // the attribute which defines the display name instead of the var's name
char *title;                // Title of the (single) plot
char *subtitle;             // subtitle of the (single) plot
char *xaxisname;            // force to plot this as X axis string
FileType xfiletype;
FileType tfiletype;
FileType vfiletype;
double xmin, xmax, vmin, vmax;
int  outidxWidth;           // how many characters should the index in output name occupy (padded with 0)
int  sizeX, sizeY;          // output image size
int  fgColor[3];            // output image foreground color rgb array (0..255 values)
int  bgColor[3];            // output image background color rgb array (0..255 values)
int  timeFix1;              // time value used for 1D 
int  timeFix2;              // time value used for 2D 
                            //   used only when there is no time iteration

// Flags from arguments or defaults
bool output_xml;
bool output_separate_dirs;
bool output_agr;
bool use_global_minmax;    // plot with calculate global minmax values for each y variable
bool skip_undef;           // skip a plot if all values are undef
bool use_regexp;           // use varmasks as regular expressions
bool list_only;            // list the variables and exit
bool multiplot;            // Make only one image that draws all variables (and all times!!!)
bool multiplot_time;       // Make one image per variable containing all timesteps
bool hide_axes;            // Do not plot the axes (2D only yet but maybe will be used for 1D)

// Variables related to 2D graphics
char *yvar;                 // name of y variable (2D plots)
char *yvar1;                // name of y variable (2D plots)
char *yfile;                // file containing y variable for 2D structured grid plots
char *ystart;               // dimension spec starting points
char *ycount;               // dimension spec counts
char *yaxisname;            // force to plot this as Y axis string
double ymin, ymax;
FileType yfiletype;
bool plot_in_2D;           // true if 2D graphics is requested for 2D arrays
#ifdef VTK
char *colormap;            // set if given as argument
vtk_options vtkopts = VTK_OPTIONS_NEW;
bool do_trianglemesh;      // the data is on a triangle mesh, x should be the vertices and y should be the triangles
bool do_rectilinear;       // plot a rectilinear grid, x and y should be 1D arrays of coordinates
bool use_polar_coordinates;// x/y is given as radius/angle or angle/radius for a polar plot
bool do_squareplots;       // scale x or y axis in 2D to get a square plot
bool transpose;            // use transpose of v in 2D plots
#endif

// Other global variables
int err;
Array x1;// = Array_new(undefArray); // we do not know the type yet
//Array y1;// = Array_new(undefArray); // we do not know the type yet
Array x2;// = Array_new(undefArray); // we do not know the type yet
Array y2;// = Array_new(undefArray); // we do not know the type yet
char *prgname; /* argv[0] */
long timeFrom1, timeSteps1, timeDimidx1, timeStepping1;
long timeFrom2, timeSteps2, timeDimidx2, timeStepping2;
int  xistart1[MAX_DIMS], xicount1[MAX_DIMS];
int  xistart2[MAX_DIMS], xicount2[MAX_DIMS];
int  yistart1[MAX_DIMS], yicount1[MAX_DIMS];
int  yistart2[MAX_DIMS], yicount2[MAX_DIMS];
int  tistart1[MAX_DIMS], ticount1[MAX_DIMS]; 
int  tistart2[MAX_DIMS], ticount2[MAX_DIMS]; 
int  istart1[MAX_DIMS], icount1[MAX_DIMS];
int  istart2[MAX_DIMS], icount2[MAX_DIMS];
int  n1Dvars = 0; // number of variables actually plotted in 1D; used for multiplot
struct timeval tp;

// Global code variables
// For 2D code
int nphi; 
int first_data_read;
Array nphi_v;
// For 1D code
// Length of psi to know which variables to plot against psi
// Plot other 1d vars vs a newly created array from 1-npoints
int psi_length;
// File to watch and file attributes (mostly file size)
char * oneddiagfile;
struct stat fileattrib;
long int previous_file_size;
// Number of steps in file and last time step processed
int tsteps;
int last_step_processed;
int current_step_processed;
int number_steps_processed;
// Global ADIOS vars
static ADIOS_GROUP * g = 0;
static ADIOS_FILE * f = 0;

//#define ESIMON
#undef ESIMON
#ifdef ESIMON
/* Include dashboard header file */
#include "esimmon.h"
#else
struct estruct { char *rundir; };
#endif
#ifdef ESIMON
ESIMMON *e;
#else
struct estruct e1 = { "." };
struct estruct * e = &e1;
#endif

// Prototypes
int find_next_step_string(char *str, char *orig);
int find_next_step(char *str, char **retstr);
static void parseDimSpec(char *str, int *dims, bool allow_time, long *tDimidx, long *tValue, long *tStepping);
static void parseRGB(char *str, int rgb[3]);

int opt;

struct option options[] = {
    {"help",                 no_argument,          NULL,    'h'},
    {"debug",                no_argument,          NULL,    'd'},
    {"input",                required_argument,    NULL,    'i'},
    {"output",               required_argument,    NULL,    'o'},
    {"shotname",             required_argument,    NULL,    's'},
    {NULL,                   0,                    NULL,    0}
};


static const char *optstring = "-do:i:s:";

void init_globals(void) {
    int i;
    // variables for arguments
    shotname             = NULL;
    outpath              = NULL;
    inpath               = NULL;
    xvar                 = NULL;
    xvar1                 = NULL;
    xaxisname            = NULL;
    tvar                 = NULL;
    xfile                = NULL;
    tfile                = NULL;
    vfile                = NULL;
    xfiletype            = undef_file;
    tfiletype            = undef_file;
    vfiletype            = undef_file;
    minmaxfile           = NULL;
    start                = NULL;
    count                = NULL;
    xstart               = NULL;
    xcount               = NULL;
    nameattr             = NULL;
    title                = NULL;
    subtitle             = NULL;
    verbose              = 0;
    list_only            = false;
    output_xml           = false;
    output_separate_dirs = false;
    output_agr           = false;
    use_global_minmax    = false;
    skip_undef           = false;
    use_regexp           = false;
    print_provenance     = false;
    multiplot            = false;
    multiplot_time       = false;
    hide_axes            = false;
    xmin                 = NAN;
    xmax                 = NAN;
    vmin                 = NAN;
    vmax                 = NAN;
    //other
    timeFrom1            = 0;
    timeFrom2            = 0;
    timeSteps1           = 0;
    timeSteps2           = 0;
    timeStepping1        = 1;
    timeStepping2        = 1;
    timeDimidx1          = -1;
    timeDimidx2          = -1;
    timeFix1             = -1;
    timeFix2             = -1;
    outidxWidth          = 4;
    sizeX                = 0;
    sizeY                = 0;
    fgColor[0]           = 0;
    fgColor[1]           = 0;
    fgColor[2]           = 0;   // default is black foreground
    bgColor[0]           = 255;
    bgColor[1]           = 255;
    bgColor[2]           = 255; // default is white background
    // set defaults
    for (i=0; i<MAX_DIMS; i++) {
        istart1[i]  = 0;
        icount1[i]  = -1;  // read full var by default
        xistart1[i] = 0;
        xicount1[i] = -1;  // read full x var by default
        tistart1[i] = 0;
        ticount1[i] = -1;  // read full t var (cannot be changed by arguments)
        istart2[i]  = 0;
        icount2[i]  = -1;  // read full var by default
        xistart2[i] = 0;
        xicount2[i] = -1;  // read full x var by default
        tistart2[i] = 0;
        ticount2[i] = -1;  // read full t var (cannot be changed by arguments)
    }

    // 2D  (variables are initialized even for not #ifdef VTK case)
    plot_in_2D          = false;  // in plotter
    yvar                = NULL;
    yvar1               = NULL;
    yaxisname           = NULL;
    yfile               = NULL;
    yfiletype           = undef_file;
    ystart              = NULL;
    ycount              = NULL;
    ymin                = NAN;
    ymax                = NAN;
    for (i=0; i<MAX_DIMS; i++) {
        yistart1[i] = 0;
        yicount1[i] = -1;  // read full y var by default
        yistart2[i] = 0;
        yicount2[i] = -1;  // read full y var by default
    }
#ifdef VTK
    plot_in_2D            = true;  // in plotter2D
    colormap              = NULL;
    do_trianglemesh       = false;
    do_rectilinear        = false;
    use_polar_coordinates = false;
    do_squareplots        = false;
    transpose             = false;
    vtk_init();
#endif
}


bool isDefined(char *var, char* name, bool printerror) {
    /* Check if required variable is defined */
    if (var == NULL) {
        if (printerror)
            fprintf(stderr, "Error: %s is not defined!\n", name);
        return false;
     }

    return true;
}

/** Check for missing args and contradictions.
  */
bool checkArgs(void) {
    bool retval = true;

    // we need an input directory file
    retval &= isDefined(inpath, "inpath", true);

    // we need an output directory file
    retval &= isDefined(outpath, "outpath", true);

    // we need an output directory file
    retval &= isDefined(shotname, "shotname", true);

    return retval;
}

/** Set default values for those variables that are not set and need not be set with an argument
  */
static void setDefaultsAfterArgs(void) {
#ifdef XMGRACE
    if (sizeY > 0) grace_setImageSize(sizeX, sizeY);
    grace_setForeground(fgColor);
    grace_setBackground(bgColor);
#endif
#ifdef VTK
    if (sizeY > 0) vtk_setImageSize(sizeX, sizeY);
    vtk_setForeground(fgColor);
    vtk_setBackground(bgColor);
#endif
}

/** Set the start[] and count[] arrays to read variables.
  * start:  "[t][-]n[,...]"; 't' in a dimension denotes the time dimension
  *                          only one t is allowed
  * count:  "[-]n[,...]"; -1 denotes 'until end'
  */
static void processDimSpecs(void) {
    // Use SCHEMA 
    //TODO: hardcoding some info based on plotter cmd for XGC1 workflow
    // will need to be changed & generalized
    timeDimidx1 = 0;
    timeDimidx2 = 1;
    // output in separate directories
    output_separate_dirs = true;
    // process var
    if (timeDimidx1 > -1) {
        timeFrom1 = istart1[ timeDimidx1 ];
        timeSteps1 = icount1[ timeDimidx1 ];
    }
    if (timeDimidx2 > -1) {
        timeFrom2 = istart2[ timeDimidx2 ];
        timeSteps2 = icount2[ timeDimidx2 ];
    }
    //if (xstart != NULL) parseDimSpec(xstart, xistart, false, NULL, NULL, NULL);
    //if (xcount != NULL) parseDimSpec(xcount, xicount, false, NULL, NULL, NULL);

#ifdef VTK
    // process y
    //if (ystart != NULL) parseDimSpec(ystart, yistart, false, NULL, NULL, NULL);
    //if (ycount != NULL) parseDimSpec(ycount, yicount, false, NULL, NULL, NULL);
#endif
}

#define PRINT_DIMS(str, v) printf("  %s = { ", str); \
    for (i=0; i<MAX_DIMS;i++) printf("%d ", v[i]);    \
    printf("}\n");

void determineMinMax( char *varname, Array a, double *minval, double *maxval) {
    double gmin, gmax;

    if (use_global_minmax) {
        // read the global min/max from file
        minmax_get(varname, Array_get_undef_value(a), &gmin, &gmax);
    }

    // determine the min value
    if (!isnan(vmin))
        *minval = vmin;        // given as argument
    else if (use_global_minmax)
        *minval = gmin;        // got value from minmax file
    else {
        *minval = Array_min(a); // calculate from the array
    }

    // determine the max value
    if (!isnan(vmax))
        *maxval = vmax;        // given as argument
    else if (use_global_minmax)
        *maxval = gmax;        // got value from minmax file
    else {
        *maxval = Array_max(a); // calculate from the array
    }

}

#ifdef VTK
/** 2D plot a 2D array with vtk.
  * Input v is expected to be a 2D array except for triangle meshes where v must be a 1D array.
  * min_val and max_val are the min-max for the colormap of the plots
  */
int plotVar2D( char *title, Array v, char* imagefilename, double min_v, double max_v,
               Array x, char *xname, double min_x, double max_x,
               Array y, char *yname, double min_y, double max_y)
{
    int i, retval=0;
    // temp vars for error checking
    int ndims_v = 0, ndims_x = 0, ndims_y = 0;
    int dimv[2] = {0,0};
    int dimx[2] = {0,0};
    int dimy[2] = {0,0};

    // determine dimensions of all arrays
    ndims_v = Array_getAllDimensions(v, dimv);
    ndims_x = Array_getAllDimensions(x, dimx);
    ndims_y = Array_getAllDimensions(y, dimy);


    /* Note:
       structured points: x and y are undef
       rectilinear grid:  x and y are 1D arrays, their size can be different
       structured grid:   ndims of y = ndims of x, and all dimensions of y = all dimensions of x
       triangle mesh:     v is 1D, len(v)=N, x is Nx2 or Nx3 array, y is Mx3 array
    */

    // check for error: x can be used only with y
    // if only x is given, this may happen for 1D plots, and should not use it here for 2D plots
    if ( Array_getType(x) != undefArray && Array_getType(y) == undefArray ) {
        fprintf(stderr, "Warning: x cannot be used without y in 2D plots. Plot without x...\n");
        x = Array_new(undefArray);
    }

    vtkopts.title = title;
    vtkopts.imagefilename = imagefilename;
    vtkopts.xname = xname;
    vtkopts.yname = yname;
    vtkopts.auto_zoom = false;
    vtkopts.show_axes = false;
    vtkopts.do_square = false;

    if (use_polar_coordinates) {
        /**************/
        /* POLAR PLOT */
        /**************/
        //  check: v must be 2D while x and y must be 1D
        if (ndims_v != 2 || ndims_x != 1 || ndims_y != 1) {
            fprintf(stderr, "Error: For polar 2D plots, v should be a 2D array while x and y should be 1D arrays\n");
            return 1;
        }
        if (vtkopts.radius_is_first_dim) {
            //  check: 1st dim of v = length of r and  2nd dim of v = length of phi
            if ( dimx[0] != dimv[0] ) {
                fprintf(stderr, "Error: size of x must = 1st dim of v. Use --polar=y if v is ordered oppositely.\n");
                return 1;
            }
            if ( dimy[0] != dimv[1] ) {
                fprintf(stderr, "Error: size of y must = 2nd dim of v. Use --polar=y if v is ordered oppositely.\n");
                return 1;
            }
        } else {
            //  check: 2nd dim of v = length of r and  1st dim of v = length of phi
            if ( dimx[0] != dimv[1] ) {
                fprintf(stderr, "Error: size of x must = 2nd dim of v. Use --polar=x if v is ordered oppositely.\n");
                return 1;
            }
            if ( dimy[0] != dimv[0] ) {
                fprintf(stderr, "Error: size of y must = 1st dim of v. Use --polar=x if v is ordered oppositely.\n");
                return 1;
            }
        }
        retval = vtk_plotPolarGrid(v, min_v, max_v, x, min_x, max_x, y, min_y, max_y, vtkopts);
    }

    else if (do_trianglemesh) {
        /*****************/
        /* TRIANGLE MESH */
        /*****************/
        // check: v is 1D, x and y are 2D
        if (ndims_v != 1 || ndims_x != 2 || ndims_y != 2) {
            fprintf(stderr, "Error: For triangle mesh plots, v should be a 1D array while x and y should be 2D arrays\n");
            return 1;
        }
        // check: length(v) = 1st dim of x
        if (dimv[0] != dimx[0]) {
            fprintf(stderr, "Error: For triangle mesh plots, size of v = size of 1st dimension of x\n");
            return 1;
        }
        // check: 2nd dim of x is 2 or 3
        if (dimx[1] != 2 && dimx[1] != 3) {
            fprintf(stderr, "Error: For triangle mesh plots, x should be an Nx2 or Nx3 array of x-y[-z] co-ordinates of points\n");
            return 1;
        }
        // check: 2nd dim of y is 3
        if (dimy[1] != 3) {
            fprintf(stderr, "Error: For triangle mesh plots, y should be an Mx3 array of point triplets (triangles of points in x)\n");
            return 1;
        }
        // check: y is an int16 or int32 array
        if (Array_getType(y) != int16Array && Array_getType(y) != int32Array) {
            fprintf(stderr, "Error: y must be an array of integers (16 or 32 bit) for triangle mesh plots.\n");
            return 1;
        }
        retval = vtk_plotTriangleMesh(v, min_v, max_v, x, min_x, max_x, y, min_y, max_y, vtkopts);
    }

    else if (do_rectilinear) {
        /********************/
        /* RECTILINEAR GRID */
        /********************/
        //  check: x and y must be 1D
        if (ndims_x != 1 || ndims_y != 1) {
            fprintf(stderr, "Error: For rectilinear plots, x and y must be 1D arrays.\n");
            return 1;
        }
        //  check: v must be 2D
        if (ndims_v != 2) {
            fprintf(stderr, "Error: For rectilinear plots, v should be a 2D array.\n");
            return 1;
        }
        //  check: 1st dim of v = length of x
        if ( dimx[0] != dimv[0] ) {
            fprintf(stderr, "Error: size of x must = 1st dim of v for rectilinear plots (x and y is 1D).\n");
            return 1;
        }
        //  check: 2nd dim of v = length of y
        if ( dimy[0] != dimv[1] ) {
            fprintf(stderr, "Error: size of y must = 2nd dim of v for rectilinear plots (x and y is 1D).\n");
            return 1;
        }
        vtkopts.auto_zoom = true;
        vtkopts.show_axes = !hide_axes;
        vtkopts.do_square = do_squareplots;
        vtkopts.transposed = transpose;
        retval = vtk_plotRectilinearGrid(v, min_v, max_v, x, min_x, max_x, y, min_y, max_y, vtkopts);
    }

    else if (ndims_v == 1 && ndims_x == 1 && ndims_y == 1) {
        /*******************/
        /* STRUCTURED GRID 1D */
        /*******************/
        //  check: dims of x = dims of v
        if (dimx[0] != dimv[0]) {
            fprintf(stderr, "Error: For structured grid 1D plots, v and x must have the same dimension.\n");
            return 1;
        }
        //  check: dims of y = dims of v
        if (dimy[0] != dimv[0]) {
            fprintf(stderr, "Error: For structured grid 1D plots, v and y must have the same dimension.\n");
            return 1;
        }
        vtkopts.auto_zoom = true;
        vtkopts.show_axes = false;
        vtkopts.do_square = do_squareplots;
        vtkopts.transposed = transpose;
        retval = vtk_plotStructuredGrid1D(v, min_v, max_v, x, min_x, max_x, y, min_y, max_y, vtkopts);
    }

    else if (ndims_x == 1 && ndims_y == 1) {
        /*********************************************************/
        /* plot a 2D array as is but with x and y as axis values */
        /*********************************************************/
        vtkopts.auto_zoom = true;
        vtkopts.show_axes = !hide_axes;
        vtkopts.do_square = do_squareplots;
        vtkopts.transposed = transpose;
        retval = vtk_plotArray(v, min_v, max_v, x, min_x, max_x, y, min_y, max_y, vtkopts);
    }

    else if (ndims_x == 0 && ndims_y == 0) {
        /*************************/
        /* plot a 2D array as is */
        /*************************/
        vtkopts.auto_zoom = true;
        vtkopts.show_axes = !hide_axes;
        vtkopts.do_square = do_squareplots;
        vtkopts.transposed = transpose;
        retval = vtk_plotArray(v, min_v, max_v, NULL, 0, dimv[0]-1, NULL, 0, dimv[1]-1, vtkopts);
    }

    else {
        fprintf(stderr,"Error: Could not figure out what kind of 2D plot to make from your data:\n");
        fprintf(stderr,"       v: [%d ", dimv[0]);
        for (i=1; i<ndims_v; i++) fprintf(stderr, ", %d", dimv[i]);
        fprintf(stderr,"]  x: [ ");
        for (i=0; i<ndims_x; i++) fprintf(stderr, " %d", dimx[i]);
        fprintf(stderr,"]  y: [ ");
        for (i=0; i<ndims_y; i++) fprintf(stderr, " %d", dimy[i]);
        fprintf(stderr,"]\n");
        fprintf(stderr,"       Options -polar or -triangle-mesh was NOT given\n");
        return 1;
    }

    return retval;
}
#endif

/* Plotting function for 2D vars */
int plot2DVariables( DataFile dfv, VarInfo vi, char* outspec,
              Array x, char *xname, double min_x_in, double max_x_in,
              Array y, char *yname, double min_y_in, double max_y_in,
              Array t, char *tvar,
              char * datafile
              )
{
    int iter, tb, te, retval=0;
    char *vname = vi.dispname;
    char *plottitle = title;
    char *legend = NULL;
    char tmp[128], outidx[32], timestr[32], outnamespec[256],newvname[256],newpath[256],planename[256];
    bool should_free_x = false;
    double min_v, max_v, min_x, max_x; //, min_y, max_y;
    Array v = Array_new(undefArray);
#ifdef VTK
    Array vT;
#endif

    if (xname && xname[0] == '/') xname++; // get rid of starting / in names
    if (yname && yname[0] == '/') yname++;
    if (vname && vname[0] == '/') vname++;

    // set iterative dimension
    // Use SCHEMA
    if (timeDimidx2 > -1 && strcmp(vname,"pot0")) {
        // negative index =  last+index
        if (timeFrom2 < 0) tb = vi.dimsize[timeDimidx2] + timeFrom2 % vi.dimsize[timeDimidx2] + dfv.time_start;
        else tb = timeFrom2;
        if (timeSteps2 < 0) te = vi.dimsize[timeDimidx2] + timeSteps2 % vi.dimsize[timeDimidx2] + dfv.time_start;
        else te = tb + timeSteps2*timeStepping2 - 1;
        // error checking
        if (tb >= vi.dimsize[timeDimidx2] + dfv.time_start) {
            fprintf(stderr,"Warning: time begin index is > steps available for variable %s\n", vi.name);
            return 0;  // not an error but cannot make plots for this variable
        }
        if (te >= vi.dimsize[timeDimidx2] + dfv.time_start) {
            fprintf(stderr,"Warning: time end index is > steps available for variable %s\n", vi.name);
            te = vi.dimsize[timeDimidx2] - 1;
        }
         if (verbose) printf("timeFrom = %ld  timeSteps = %ld  timeStepping = %ld -> from = %d  to = %d\n",
             timeFrom2, timeSteps2, timeStepping2, tb, te);
    } else {
        // no iter dimension was specified
        tb = 0;
        te = 0;
    }

    if (verbose) printf("variable %s ...\n", vname);
    for (iter = tb; iter <= te; iter+=timeStepping2) {
        if (verbose) printf("  iteration %d ...\n",iter);
        if (timeDimidx2 > -1 && strcmp(vname,"pot0")) {
            istart2[timeDimidx2] = iter;
            icount2[timeDimidx2] = 1;
            // make the output name idx string based on iter and variable t
            if (t != NULL && Array_getType(t) != undefArray) {
                Array_strvalue(t, iter, timestr);
                Array_strvalue_padded(t, iter, outidx, outidxWidth);
                if (verbose>1) printf("    Output file index: %s\n", timestr);
            } else {
                // So here we have a case of tricking the time concept
                // timeFix2 is used for time step of the simulation and
                // file output like xgc.3d.timestep.bp
                // But for each variable we have 32 places in each file
                // that is not a time dimension but we are going to use it
                // as so to create diffrerent variable names
                snprintf(timestr, 32, "%d", iter);
                snprintf(outidx, 32, "%0*d", outidxWidth, iter);
            }
            snprintf(timestr, 32, "%0*d", outidxWidth, timeFix2);
            snprintf(planename,256, "%s-%s", outspec,outidx);
            snprintf(newvname,256, "%s/%s-%s", vname,vname,outidx);
            snprintf(newpath,256, "%s/%s/%s", outpath,vname,newvname);
            createdir_recursive(planename);
            snprintf(outnamespec, 256, "%s/%s-%s.%s", planename, vname, outidx, timestr);
        } else if (timeFix2 > -1) { // timefix counts only if there is no iteration
            snprintf(timestr, 32, "%d", timeFix2);
            snprintf(newvname,256, "%s",vname);
            snprintf(outnamespec, 256, "%s.%0*d", outspec, outidxWidth, timeFix2);
        } else {
            strncpy(outnamespec, outspec, 256);
        }

#ifdef ESIMON 
            err = esimmon_addVariable(e, newvname, timeFix2,
                    datafile, esimmon_transfer_external, NULL, 0, "");
            if (err) {
                fprintf(stderr, "ERROR at addVariable: errno=%d\n%s\n",
                        esimmon_errno, esimmon_get_last_errmsg());
            }
#endif

        _tstart();
        v = reader_read_var(dfv, vi, istart2, icount2);
        _tend();
        _tdiff(&tp);
        if (verbose>1) printf("     time to read %s = %d.%6.6d sec\n", vi.name, tp.tv_sec, tp.tv_usec);
        if ( Array_errno(v) ) {
            fprintf(stderr, "Error: could not read var %s, iteration %d from %s\n", vi.name, iter, vfile);
            exit(5);
        }
        if (verbose) printf("    v array length = %d\n", Array_length(v));
#ifdef VTK
        // if --transpose then transpose v now
        if (transpose && Array_getDimensions(v) == 2) {
            if (verbose) printf("   Transpose variable %s ...\n", vname);
            vT = Array_transpose(v);
            if ( Array_getType(vT) == undefArray) {
                return 2;
            }
            Array_free(v);
            v = vT;
        }
#endif

#ifdef VTK
        determineMinMax(vname, v, &min_v, &max_v);
        if (plottitle == NULL) {
            snprintf(tmp, 128, "%s", vname);
            plottitle = tmp;
        }
        retval = plotVar2D( plottitle, v, outnamespec, min_v, max_v,
                x, xname, min_x_in, max_x_in,
                y, yname, min_y_in, max_y_in);
        // write minmax into file
        if (!retval && minmaxfile != NULL && !use_global_minmax) {
            if (verbose>1) printf("    print minmax\n");
            minmax_print(vname, iter, Array_min(v), Array_max(v));
        }
        Array_free(v);
        if (retval) break;  // exit after first error
        continue;
#endif

        // if x is empty, let's create a [0..n-1] array here
        if (x == NULL || Array_getType(x) == undefArray) {
            min_x = 0;
            max_x = Array_length(v)-1;
            x = Array_indgen((int)min_x,(int)max_x,1);
            should_free_x = true;
            if (verbose>1) printf("    created x array from [%d..%d], len=%ld, type=%s\n",
                (int)min_x, (int)max_x, Array_length(x), arrayTypeStrings[Array_getType(x)]);
            // set to min/max if provided as argument
            if (!isnan(min_x_in)) min_x = min_x_in;
            if (!isnan(max_x_in)) max_x = max_x_in;
            if (verbose>1) printf("    x min/max %d..%d\n", (int)min_x, (int)max_x);
        } else {
            min_x = min_x_in;
            max_x = max_x_in;
        }

        // write minmax into file
        if (minmaxfile != NULL && !use_global_minmax) {
            if (verbose>1) printf("    print minmax\n");
            minmax_print(vname, iter, Array_min(v), Array_max(v));
        }

        // free array(s) up for the next iteration
        Array_free(v);
        if (should_free_x) {
            Array_free(x);
            should_free_x = false;
            x = NULL;
        }
    }

    if (verbose>1) printf("    return from plot2DVariables\n");
    return retval;
}

/* Plot 1D variables */
int plot1DVariables( DataFile dfv, VarInfo vi, char* outspec,
              Array x, char *xname, double min_x_in, double max_x_in,
              Array y, char *yname, double min_y_in, double max_y_in,
              Array t, char *tvar
              )
{
    int iter, tb, te, retval=0;
    char *vname = vi.dispname;
    char *plottitle = title;
    char *legend = NULL;
    char tmp[128], outidx[32], timestr[32], outnamespec[256];
    bool should_free_x = false;
    double min_v, max_v, min_x, max_x; //, min_y, max_y;
    Array v = Array_new(undefArray);
#ifdef VTK
    Array vT;
#endif

    if (xname && xname[0] == '/') xname++; // get rid of starting / in names
    if (yname && yname[0] == '/') yname++;
    if (vname && vname[0] == '/') vname++;

    // set iterative dimension
    if (timeDimidx1 > -1) {
        // negative index =  last+index
        if (timeFrom1 < 0) tb = vi.dimsize[timeDimidx1] + timeFrom1 % vi.dimsize[timeDimidx1] + dfv.time_start;
        else tb = timeFrom1;
        if (timeSteps1 < 0) te = vi.dimsize[timeDimidx1] + timeSteps1 % vi.dimsize[timeDimidx1] + dfv.time_start;
        else te = tb + timeSteps1*timeStepping1 - 1;
        // error checking
        if (tb >= vi.dimsize[timeDimidx1] + dfv.time_start) {
            fprintf(stderr,"Warning: time begin index is > steps available for variable %s\n", vi.name);            return 0;  // not an error but cannot make plots for this variable
        }
        if (te >= vi.dimsize[timeDimidx1] + dfv.time_start) {
            fprintf(stderr,"Warning: time end index is > steps available for variable %s\n", vi.name);
            te = vi.dimsize[timeDimidx1] - 1;
        }
         if (verbose) printf("timeFrom = %ld  timeSteps = %ld  timeStepping = %ld -> from = %d  to = %d\n",
             timeFrom1, timeSteps1, timeStepping1, tb, te);
    } else {
        // no iter dimension was specified
        tb = 0;
        te = 0;
    }
    if (verbose) printf("variable %s ...\n", vname);
    tb = last_step_processed;
    number_steps_processed = 0;
    for (iter = tb; iter <= te; iter+=timeStepping1) {
        if (verbose) printf("  iteration %d ...\n",iter);
        if (timeDimidx1 > -1) {
            istart1[timeDimidx1] = iter;
            icount1[timeDimidx1] = 1;
            number_steps_processed++;
            // make the output name idx string based on iter and variable t
            if (t != NULL && Array_getType(t) != undefArray) {
                Array_strvalue(t, iter, timestr);
                Array_strvalue_padded(t, iter, outidx, outidxWidth);
                if (verbose>1) printf("    Output file index: %s\n", timestr);
            } else {
                snprintf(timestr, 32, "%d", iter);
                snprintf(outidx, 32, "%0*d", outidxWidth, iter);
            }
            snprintf(outnamespec, 256, "%s.%s", outspec, outidx);
        } else if (timeFix1 > -1) { // timefix counts only if there is no iteration
            snprintf(timestr, 32, "%d", timeFix1);
            snprintf(outnamespec, 256, "%s.%0*d", outspec, outidxWidth, timeFix1);
        } else {
            strncpy(outnamespec, outspec, 256);
        }
#ifdef ESIMON 
            // Add to ESIMON 
            err = esimmon_addVariable(e, vname, iter,
                    oneddiagfile, esimmon_transfer_external, NULL, 0, "");
            if (err) {
                fprintf(stderr, "ERROR at addVariable: errno=%d\n%s\n",
                        esimmon_errno, esimmon_get_last_errmsg());
            }
#endif
        // read in variable (Array v will be allocated)
        _tstart();
        v = reader_read_var(dfv, vi, istart1, icount1);
        _tend();
        _tdiff(&tp);
        if (verbose>1) printf("     time to read %s = %d.%6.6d sec\n", vi.name, tp.tv_sec, tp.tv_usec);
        if ( Array_errno(v) ) {
            fprintf(stderr, "Error: could not read var %s, iteration %d from %s\n", vi.name, iter, vfile);
            exit(5);
        }
        if (verbose) printf("    v array length = %d\n", Array_length(v));

        // if x is empty, let's create a [0..n-1] array here 
        if (x == NULL || Array_getType(x) == undefArray) {
            min_x = 0;
            max_x = Array_length(v)-1;
            x = Array_indgen((int)min_x,(int)max_x,1);
            should_free_x = true;
            if (verbose>1) printf("    created x array from [%d..%d], len=%ld, type=%s\n",                (int)min_x, (int)max_x, Array_length(x), arrayTypeStrings[Array_getType(x)]);
            // set to min/max if provided as argument
            if (!isnan(min_x_in)) min_x = min_x_in;
            if (!isnan(max_x_in)) max_x = max_x_in;
            if (verbose>1) printf("    x min/max %d..%d\n", (int)min_x, (int)max_x);
        } else {
            min_x = min_x_in;
            max_x = max_x_in;
        }
#ifdef XMGRACE        
            // make plot
            if (iter==tb) n1Dvars++; // this variable is plotted in 1D
            // .png, .agr will be added in grace_plot to outnamespec
            determineMinMax(vname, v, &min_v, &max_v);
            if (!skip_undef || !Array_is_data_undef_all(v))
            {
                if (plottitle == NULL) {
                    if (xname) snprintf(tmp, 128, "%s vs. %s", vname, xname);
                    else snprintf(tmp, 128, "%s", vname);
                    if (verbose) printf("    plot %s  min=%le  max=%le\n", outnamespec, min_v, max_v);
                    plottitle = tmp;
                }

                // init plot now if needed (3 cases: multiplot, multiplot_time, normal) 
                if (multiplot) {
                    // error check (should not reach because checkArgs() should have stopped already)
                    if (iter != tb) {
                        fprintf(stderr, "Error: Plotter cannot plot several variables into one image _and_ loop over time to generate many images.\n");
                        fprintf(stderr, "In case of --multiplot, do not use time spec in --start.\n");                        fprintf(stderr, "Use --multiplot-time to plot several timesteps of one variable on one plot.\n");
                        exit(180);
                    }
                    // init for first var and first time only
                    if (n1Dvars == 1 && iter==tb) {
                        if (title == NULL) {
                            snprintf(tmp, 128, "Multiplot");
                            plottitle = tmp;
                        }
                        grace_init_plot (plottitle, subtitle, xname, vname, 0.30);
                    }
                    legend = vname;
                } else if (multiplot_time) {
                    // init plot at the first timestep
                    if (iter==tb) {
                        grace_init_plot (plottitle, subtitle, xname, vname, 0.157051);
                    }
                    legend = timestr;
                } else {
                    grace_init_plot (plottitle, subtitle, xname, vname, 0.0);
                    legend = NULL;
                }

                // draw this variable and time
                grace_draw_plot( Array_length(v), x, v, vname, timestr, legend, min_x, max_x, min_v, max_v );

                // save image if not multiplot_time or multiplot
                if (!multiplot_time && !multiplot){
                    grace_save_plot(outnamespec, output_agr);
                }
            } else {
                fprintf(stderr, "    plot %s is skipped because values are undefined\n", outnamespec);
            }
#else /* XMGRACE not available */
            fprintf(stderr, "ERROR: this code is not compiled with xmgrace, so it cannot make 1D plots. Either do 2D plots, or print text with -p or --output-xml\n");
            break; /* break the time loop */

#endif

        // write minmax into file
        if (minmaxfile != NULL && !use_global_minmax) {
            if (verbose>1) printf("    print minmax\n");
            minmax_print(vname, iter, Array_min(v), Array_max(v));
        }

        // free array(s) up for the next iteration
        Array_free(v);
        if (should_free_x) {
            Array_free(x);
            should_free_x = false;
            x = NULL;
        }
    }

#ifdef XMGRACE
    // save image if multiplot_time 
    if (multiplot_time) {
        // use outspec to save instead of outnamespec which has time in it
        grace_save_plot(outspec, output_agr);
    }
#endif

    if (verbose>1) printf("    return from plotVar\n");
    return retval;
}

static bool conformsToDimSpec(VarInfo vi) {
    if (timeDimidx2 > vi.ndims) {
            fprintf(stderr,"The time dimension (%ld) is out of the dimensions of this variable %s (%d). Skip it.\n",
                    timeDimidx2, vi.name, vi.ndims);
            return false;
    }
    if (timeDimidx2 == 0 && vi.ndims == 1) {
            fprintf(stderr,"Variable %s has only 1 dimension and time is specified. Skip it.\n", vi.name);
            return false;
    }
    return true;
}

/* Check file info before doing work */
int checkFile(){
    stat(oneddiagfile, &fileattrib);
    ushort mode = fileattrib.st_mode;
    // printf("Compare: %ld and %ld\n", previous_file_size, fileattrib.st_size);
    // FIXME if (fileattrib.st_size == previous_file_size){
    //    printf("file has not increased, wait\n");
    //    return 0;
    //}else{
        // printf("Check time step %d\n", last_step_processed);
        // Check ADIOS file
        f = adios_fopen (oneddiagfile, MPI_COMM_WORLD);
        if (f == NULL) {
            fprintf (stderr, "%s\n", adios_errmsg());
            return 1;
        }

        // Open the group - don't need the group name, can use grouname[0]
        g = adios_gopen (f, f->group_namelist[0]);
        if (g == NULL) {
            fprintf (stderr, "%s\n", adios_errmsg());
            return 1;
        }

        // Check time variable in adios file
        // SPECIFIC to XGC1 - USE SCHEMA
        ADIOS_VARINFO *v0 = adios_inq_var (g, "time");
        // printf("length of this var: %d\n",v0->ndim);
        // printf("dim1: %d\n",v0->dims[0]);
        // printf("now compare %d and %d\n", last_step_processed, v0->dims[0]);
        if (last_step_processed < v0->dims[0]){
            current_step_processed = v0->dims[0];
            printf("Process all new time steps until t=%d\n", current_step_processed);
            do1DWork();
        }else{
            printf("Nothing new to process.\n");
        }
        adios_free_varinfo(v0);
        adios_gclose (g);
        adios_fclose (f);
    //}
    return 0;

}

/** Main body of work for 1D vars. **/
int do1DWork() {
    char buf[256], *varbasename; // to get base name of a variable with full path
    DataFile df_v = DataFile_new();
    DataFile df_x = DataFile_new();
    DataFile df_y = DataFile_new();
    DataFile df_t = DataFile_new();
    Array x = Array_new(undefArray); // we do not know the type yet
    Array y = Array_new(undefArray); // we do not know the type yet
    Array t = Array_new(undefArray); // we do not know the type yet
    Array time = Array_new(undefArray); // we do not know the type yet
    char outspec[256];
    int i, excode;
    int ntsteps;
    int nVarsMatched = 0;
    double min_x = xmin, max_x = xmax, min_y = ymin, max_y = ymax;
#ifdef VTK
    int ndims_x, ndims_y, dimx, dimy; // temp vars for error checking
#endif
    bool outpath_is_dir;
    int retval = 0;

    outpath_is_dir  = is_dir(outpath);  // so we run stat only once (+ in checkArgs)

    // open yfile
    _tstart();
    df_v = reader_open(oneddiagfile, adiosbp_file, nameattr);
    _tend();
    _tdiff(&tp);
    if (verbose>1) printf("     time to open %s = %d.%6.6d sec\n", vfile, tp.tv_sec, tp.tv_usec);
    if (df_v.fd == -1 && df_v.fh == -1) {
        fprintf(stderr, "Error: could not open file %s\n", vfile);
        retval=1; goto finalize;
    }

    // xvar is psi
    // SPECIFIC to XGC1 - USE SCHEMA
    xvar = strdup("/psi");
    // FIXME psi does not really change but it could
    // xstart should move with the time step
    // Schema will fix guessing
    xicount1[0] = 1;
    xicount1[1] = -1;

    // plot any var which second dimension is psi against psi
    // plot the rest using x=null to just plot the 1d array
    // get the time array to check time steps
    // read x var if should be read from xfile
    _tstart();
    df_x = reader_open(oneddiagfile, adiosbp_file, nameattr);
    _tend();
    _tdiff(&tp);
    if (verbose>1) printf("     time to open %s = %d.%6.6d sec\n", xfile, tp.tv_sec, tp.tv_usec);
    if (df_x.fd == -1 && df_x.fh == -1) {
        fprintf(stderr, "\nError: could not open file %s for x var %s\n", xfile, xvar);
        retval=2; goto finalize;
    }

    // read in x variable (Array x will be allocated)
    _tstart();
    x = reader_read_var_byname(df_x, xvar, xistart1, xicount1);
    _tend();
    _tdiff(&tp);
    if (verbose>1) printf("     time to read x %s = %d.%6.6d sec\n", xvar, tp.tv_sec, tp.tv_usec);
    if ( Array_errno(x) ) {
        fprintf(stderr, "\nError: could not read x var %s from %s\n", xvar, xfile);
        retval=3; goto finalize;
    }
    if (isnan(min_x)) min_x = Array_min(x); // if xmin is not argument, get from array
    if (isnan(max_x)) max_x = Array_max(x);

    psi_length = Array_length(x);

    // read t var if should be read from tfile
    if (tvar != NULL) {
        // open tfile, if not the same as vfile
        if (strcmp(tfile, vfile)) {
            df_t = reader_open(tfile, tfiletype, nameattr);
            if (df_t.fd == -1 && df_t.fh == -1) {
                fprintf(stderr, "\nError: could not open file %s for t var %s\n", tfile, tvar);
                retval=8; goto finalize;
            }
        } else {  // cannot be both NULL, because we check that in checkArgs
            if (verbose>2) printf("t file is same as var file %s\n", vfile);
            df_t = df_v;
        }

        // read in t variable (Array t will be allocated)
        t = reader_read_var_byname(df_t, tvar, tistart1, ticount1);
        if ( Array_errno(t) ) {
            fprintf(stderr, "\nError: could not read t var %s from %s\n", tvar, tfile);
            retval=9; goto finalize;
        }
    } // t remains undefArray if not read in here

    // Open minmax file, so that codes deeper in the calls can write/read minmax data
    if (minmaxfile)
        minmax_open(minmaxfile, (use_global_minmax ? 1 : 0));

    // Determine how many vars are in this file
    // Determine the last time step processed
    // Determine the time step to process
    // printf("There are: %d variables in this file.\n", df_v.nvars);
    // SPECIFIC to XGC1 - USE SCHEMA
    for (i=0; i < df_v.nvars; i++) {
        // Don't plot psi and psi_msk
        if (!strcmp(df_v.varinfo[i].name,"/psi") || !strcmp(df_v.varinfo[i].name,"/psi_mks"))
            continue;

        // Don't plot variables that are not 2D
        if (df_v.varinfo[i].ndims != 2)
            continue;

       ntsteps = df_v.varinfo[i].dimsize[0];

        excode = 0;
        nVarsMatched++;
        // get base name of the variable (can be a full path in HDF5/ADIOS)
        strncpy(buf, df_v.varinfo[i].dispname, 256);
        varbasename = basename(buf);

        // prepare output name (without numbering and extension)
        if (multiplot) {
            snprintf(outspec, 256, "%s", outpath);
        } else if (output_separate_dirs) {
            snprintf(outspec, 256, "%s/%s", outpath, df_v.varinfo[i].dispname);
            excode = createdir_recursive(outspec);
            if (!excode) {
                if (df_v.varinfo[i].dispname[0] == '/'){
                    if (outpath[strlen(outpath)-1] == '/')
                        snprintf(outspec, 256, "%s%s/%s", outpath, varbasename, varbasename);
                    else
                        snprintf(outspec, 256, "%s%s/%s", outpath, df_v.varinfo[i].dispname, varbasename);
                }else{
                    snprintf(outspec, 256, "%s/%s/%s", outpath, df_v.varinfo[i].dispname, varbasename);
                }
            }
        } else if (outpath_is_dir) {
            if (varbasename[0] == '/'){
                snprintf(outspec, 256, "%s%s", outpath, varbasename);
            }else{
                if (outpath[strlen(outpath)-1] == '/')
                    snprintf(outspec, 256, "%s%s", outpath, varbasename);
                else
                    snprintf(outspec, 256, "%s/%s", outpath, varbasename);
            }
        } else {
            if (nVarsMatched == 1)
                snprintf(outspec, 256, "%s", outpath);
            else {
                fprintf(stderr, "\nError: Plotting several variables into one output file (%s) is not supported."
                    "Specify a directory as --output \n", outpath);
                break;
            }
        }
        if (!excode) {
            if (df_v.varinfo[i].dimsize[1] == psi_length){
                retval = plot1DVariables( df_v, df_v.varinfo[i], outspec,
                        x, (xaxisname != NULL ? xaxisname : xvar), min_x, max_x,
                        y, (yaxisname != NULL ? yaxisname : yvar), min_y, max_y,
                        t, tvar);
            }else{
                retval = plot1DVariables( df_v, df_v.varinfo[i], outspec,
                        NULL, (xaxisname != NULL ? xaxisname : xvar), min_x, max_x,
                        y, (yaxisname != NULL ? yaxisname : yvar), min_y, max_y,
                        t, tvar);
            }
        } else {
            retval = excode;
        }
        if (retval) break;  // break after first error
    }
#ifdef XMGRACE
    if (multiplot && n1Dvars > 0) {
        // in case of multiplot we save the plot image here
        grace_save_plot(outpath, output_agr);
    }
#endif

    if (!retval && nVarsMatched == 0) {
        fprintf(stderr, "\nError: None of the variables matched any yvar name/regexp you provided\n");
        retval = 10;
    }

finalize:    // we jump here in case of errors too

    if (minmaxfile)    minmax_close();

    // close data file(s)
    if (tfile != NULL && strcmp(tfile, xfile) && strcmp(tfile, vfile))
        reader_close(df_t);
    if (yfile != NULL && strcmp(yfile, xfile) && strcmp(yfile, vfile))
        reader_close(df_y);
    if (xfile != NULL && strcmp(xfile, vfile))
        reader_close(df_x);
    reader_close(df_v);

    return retval;
}

/*
 * Instead of doing all the work for one variable
 * or reading the mesh for each variable, store
 * the mesh elements.
 * 
 */

int storeMesh(char *meshfile) {
    DataFile df_x = DataFile_new();
    int i; 
    double min_x = xmin, max_x = xmax; // min_y = ymin, max_y = ymax;
#ifdef VTK
    int ndims_x, ndims_y, dimx, dimy; // temp vars for error checking
#endif
    bool outpath_is_dir;
    int retval = 0;

    outpath_is_dir  = is_dir(outpath);  // so we run stat only once (+ in checkArgs)

    // read x var if should be read from xfile
    // if (xvar != NULL) {
    // open xfile, if not the same as vfile
    // printf("file type is: %d\n", adiosbp_file);
    _tstart();
    df_x = reader_open(meshfile, 3, nameattr);
    _tend();
    _tdiff(&tp);

    // Use SCHEMA    
    // read in x variable (Array x will be allocated)
    xvar = strdup("/coordinates/values");
    _tstart();
    x2 = reader_read_var_byname(df_x, xvar, xistart2, xicount2);
    _tend();
    _tdiff(&tp);
    if (verbose>1) printf("     time to read x %s = %d.%6.6d sec\n", xvar, tp.tv_sec, tp.tv_usec);
    if ( Array_errno(x2) ) {
        fprintf(stderr, "\nError: could not read x var %s from %s\n", xvar, meshfile);
        retval=3; goto finalize;
    }
    if (isnan(min_x)) min_x = Array_min(x2); // if xmin is not argument, get from array
    if (isnan(max_x)) max_x = Array_max(x2);

#ifdef VTK
    // open yfile, if not the same as xfile or vfile
    // read in y variable (Array y will be allocated)
    yvar = strdup("/cell_set[0]/node_connect_list");
    _tstart();
    y2 = reader_read_var_byname(df_x, yvar, yistart2, yicount2);
    _tend();
    _tdiff(&tp);
    if (verbose>1) printf("     time to read y %s = %d.%6.6d sec\n", yvar, tp.tv_sec, tp.tv_usec);
    if ( Array_errno(y2) ) {
        fprintf(stderr, "\nError: could not read y var %s from %s\n", yvar, meshfile);
        retval=5; goto finalize;
    }

    // must be the same number of dimensions as x (and x exist if we are at this point)
    ndims_x = Array_getDimensions(x2);
    ndims_y = Array_getDimensions(y2);
    if (ndims_x != ndims_y) {
        fprintf(stderr, "Error: # of dimensions of x (=%d) and y (=%d) must be the same for 2D plots\n",
                ndims_x, ndims_y);
        retval=6; goto finalize;
    }

    do_trianglemesh = 1;
    // if > 1 dims and not a triangle mesh is requested then
    //   the dimensions themselves should also be equal (for 1D, x and y is probably different)
    if (!do_trianglemesh && ndims_x > 1) {
        for (i=0; i<ndims_x; i++) {
            dimx = Array_getDimensionSize(x2,i);
            dimy = Array_getDimensionSize(y2,i);
            if (dimx != dimy) {
                fprintf(stderr, "Error: dimension %d of x (=%d) and y (=%d) must be the same for 2D plots\n",
                        i, dimx, dimy);
                retval=7; goto finalize;
            }
        }
    }
    // for triangle mesh, x and y should have 2 dimensions
    if (do_trianglemesh && ndims_x != 2) {
        fprintf(stderr, "Error:  x (=%d) and y (=%d) must be 2D arrays for triangle mesh plots\n",
                ndims_x, ndims_y);
        retval=6; goto finalize;
    }
#endif

    // Open minmax file, so that codes deeper in the calls can write/read minmax data
    if (minmaxfile)
        minmax_open(minmaxfile, (use_global_minmax ? 1 : 0));

finalize:    // we jump here in case of errors too

    if (minmaxfile)    minmax_close();

    // close data file(s)
    reader_close(df_x);

    return retval;

}

/* Main body of work.
 * open value files
 * plot all variables 
 */
int plotAll2DVarsForTimeStep( char* datafile)
{
    char buf[256], *varbasename; // to get base name of a variable with full path
    DataFile df_v = DataFile_new();
    DataFile df_x = DataFile_new();
    DataFile df_y = DataFile_new();
    DataFile df_t = DataFile_new();
    Array t = Array_new(undefArray); // we do not know the type yet
    char outspec[256];
    int i, excode;
    int nVarsMatched = 0;
    double min_x = xmin, max_x = xmax, min_y = ymin, max_y = ymax;
    bool outpath_is_dir;
    int retval = 0;

    outpath_is_dir  = is_dir(outpath);  // so we run stat only once (+ in checkArgs)

    // open yfile
    _tstart();
    df_v = reader_open(datafile, adiosbp_file, nameattr);
    _tend();
    _tdiff(&tp);
    if (verbose>1) printf("     time to open %s = %d.%6.6d sec\n", vfile, tp.tv_sec, tp.tv_usec);
    if (df_v.fd == -1 && df_v.fh == -1) {
        fprintf(stderr, "Error: could not open file %s\n", vfile);
        retval=1; goto finalize;
    }

    // Open minmax file, so that codes deeper in the calls can write/read minmax data
    if (minmaxfile)
        minmax_open(minmaxfile, (use_global_minmax ? 1 : 0));

    // Read nphi first, to make sure we plot all 2D vars
    if (first_data_read == 0){
        nphi_v = reader_read_var_byname(df_v, "/nphi", 0, -1);
        nphi = (int)(Array_doublevalue(nphi_v,0));
        first_data_read = 1;
    }

    // DEBUG printf("nphi: %d\n", (int)(Array_doublevalue(nphi_v,0)));
    // DEBUG printf("nphi: %d\n", (int)(Array_doublevalue(nphi_v,0)));
    // DEBUG Array_print(nphi_v);

    for (i=0; i < df_v.nvars; i++) {
        // Use SCHEMA
        //   Plotting instead any var with second dim=nphi
        //   Still tricking the code for pot0 until asking more question
        //   about the code or adding the schema
        // Good if (strcmp(df_v.varinfo[i].name,"pot0") && (df_v.varinfo[i].ndims != 2 && df_v.varinfo[i].dimsize[1] != nphi))
        if (df_v.varinfo[i].ndims == 0)
            continue;
          
        if (!conformsToDimSpec(df_v.varinfo[i]))
            continue;
        excode = 0;
        nVarsMatched++;
        // get base name of the variable (can be a full path in HDF5/ADIOS)
        strncpy(buf, df_v.varinfo[i].dispname, 256);
        varbasename = basename(buf);
        
        // prepare output name (without numbering and extension)
        if (multiplot) {
            snprintf(outspec, 256, "%s", outpath);
        } else if (output_separate_dirs) {
            snprintf(outspec, 256, "%s/%s", outpath, df_v.varinfo[i].dispname);
            excode = createdir_recursive(outspec);
            if (!excode) {
                if (df_v.varinfo[i].dispname[0] == '/'){
                    if (outpath[strlen(outpath)-1] == '/')
                        snprintf(outspec, 256, "%s%s/%s", outpath, varbasename, varbasename);
                    else
                        snprintf(outspec, 256, "%s%s/%s", outpath, df_v.varinfo[i].dispname, varbasename);
                }else{
                    snprintf(outspec, 256, "%s/%s/%s", outpath, df_v.varinfo[i].dispname, varbasename);
                }
            }
        } else if (outpath_is_dir) {
            if (varbasename[0] == '/'){
                snprintf(outspec, 256, "%s%s", outpath, varbasename);
            }else{
                if (outpath[strlen(outpath)-1] == '/')
                    snprintf(outspec, 256, "%s%s", outpath, varbasename);
                else
                    snprintf(outspec, 256, "%s/%s", outpath, varbasename);
            }
        } else {
            if (nVarsMatched == 1)
                snprintf(outspec, 256, "%s", outpath);
            else {
                fprintf(stderr, "\nError: Plotting several variables into one output file (%s) is not supported."
                    "Specify a directory as --output \n", outpath);
                break;
            }
        }
        
        if (!excode) {
            retval = plot2DVariables( df_v, df_v.varinfo[i], outspec,
                              x2, (xaxisname != NULL ? xaxisname : xvar), min_x, max_x,
                              y2, (yaxisname != NULL ? yaxisname : yvar), min_y, max_y,
                              t, tvar,datafile);

        } else {
            retval = excode;
        }
        if (retval) break;  // break after first error
    }
#ifdef XMGRACE
    if (multiplot && n1Dvars > 0) {
        // in case of multiplot we save the plot image here
        grace_save_plot(outpath, output_agr);
    }
#endif

    if (!retval && nVarsMatched == 0) {
        fprintf(stderr, "\nError: None of the variables matched any yvar name/regexp you provided\n");
        retval = 10;
    }
         

finalize:    // we jump here in case of errors too

    if (minmaxfile)    minmax_close();

    // close data file(s)
    if (tfile != NULL && strcmp(tfile, xfile) && strcmp(tfile, vfile))
        reader_close(df_t);
    if (yfile != NULL && strcmp(yfile, xfile) && strcmp(yfile, vfile))
        reader_close(df_y);
    if (xfile != NULL && strcmp(xfile, vfile))
        reader_close(df_x);
    reader_close(df_v);

    return retval;
}


// int main (int argc, char ** argv)
int main (int argc, char * argv[])
{
        err = 0;

        nphi = 0; 
        first_data_read = 0;
        nphi_v = Array_new(undefArray);

        // What needs to be passed as arguments versus what is reasonably
        // constant at least for a while in EPSi
        // The machine that is the closest related to eSiMon (ewok-web2)
        // is sith - When we add the method of directly uploading to the
        // server, the machine will need to become an argument as well
        char *machine = "sith";
        // The shot name, needs to be an argument (REQUIRED)
        // char *shotname ="4shot";
        // Reasonable that this is for XGC1 and will be later in xml
        // or config file. Leave hard coded for now.
        char *code = "xgc1";
        // Shot directory should be given on titan (or later hopper)
        // and created on sith under project directory
        // That is the input parameter
        // char * simdir;
        // Shot dir on sith needs to be constructed from project dir
        // user name and shot number
        char *shotdir; // = "/ccs/proj/e2e/rbarreto/xgc_plotter/output4";
        // These two arrays are used to store the mesh. This is definitely
        // XGC specific because we have a 2D triangular mesh but can stay
        // hardcoded as long as the dashboard is mostly used to display
        // 1D & 2D graphics.
        x1 = Array_new(undefArray); // we do not know the type yet
        //y1 = Array_new(undefArray); // we do not know the type yet
        x2 = Array_new(undefArray); // we do not know the type yet
        y2 = Array_new(undefArray); // we do not know the type yet
        // Indicate to use eSiMon dashboard or not - could be given as 
        // option as well but is fair to assume yes for EPSi
        int use_esimmon = 0; // no

        // Therefore before getting into reading the schema and figure out
        // what to plot, we need at least
        // 1. dir on titan
        // 2. project dir on sith (visible from ewok-web2)
        // 3. shot name

        // Init global variables
        init_globals();

#ifdef VTK
   error_check_sizeofbool(sizeof(char));
#endif

   prgname = strdup(argv[0]);

   /* other variables */
   char c, last_c='_';
   // int last_opt = -1;
   /* Process the arguments */
   while ((c = getopt_long_only(argc, argv, optstring, options, NULL)) != -1) {
        switch (c) {
                case 'd':
                    verbose++;
                    break;
                case 'i':
                    inpath = strndup(optarg,256);
                    break;
                case 'o':
                    outpath = strndup(optarg,256);
                    break;
                case 's':
                    shotname = strndup(optarg,256);
                    break;
                default:
                    printf("Processing default: %c\n", c);
                    break;
        } /* end switch */
        last_c = c;
   } /* end while */

   // Usage instructions
   if (argc < 6) {
       printf("Usage: %s -i <sim-input-dir> -o <proj-output-dir> -s <shotname>\n", argv[0]);
       return 1;
   }

   setDefaultsAfterArgs();
   processDimSpecs();
   if(!checkArgs())
       return 1;

   // At this point compose the output directory
   // Get username
   char *username;
   username = getenv("USER");
   shotdir = malloc(strlen(outpath)+strlen(username)+strlen(shotname)+3);
   strcpy(shotdir,outpath);
   if (outpath[strlen(outpath)-1] != '/')
       strcat(shotdir,"/");
   strcat(shotdir,username);
   strcat(shotdir,"/");
   strcat(shotdir,shotname);
   int cdir = createdir_recursive(shotdir);
   if (cdir != 0)
       printf("Error creating shot directory in project space: %s\n", outpath);

   if (use_esimmon == 0){
       outpath = strdup(shotdir);
   }else{
       outpath = strdup(e->rundir);
   }

   setDefaultsAfterArgs();
   processDimSpecs();
   if(!checkArgs())
       return 1;

#ifdef ESIMON
    char * p = getenv("ESIMMON_METHOD");
    if (p) {
        use_esimmon = 1;
    } else {
        fprintf(stderr, "WARNING: ESIMMON_METHOD env var not set. Will not use eSiMon\n");
        use_esimmon = 0;
    }
    if (use_esimmon) {
        esimmon_init(NULL, NULL);
        e = esimmon_initRun( shotname, RUN_SEPARATE, code, machine, shotdir);
        if (!e) {
            fprintf(stderr, "ERROR at initRun: errno=%d\n%s\n",
                    esimmon_errno, esimmon_get_last_errmsg());
            exit(1);
        }
        printf("Run dir=%s\n", e->rundir);
        printf("Run name=%s\n", e->runname);
        printf("Cheatsheet ID=%ld\n", e->id);
        err = createdir_recursive(e->rundir);
        if (err) {
            fprintf(stderr, "ERROR at creating directory %s, mkdir error code=%d\n\n",
                    e->rundir, err);
            exit(1);
        }
    } else {
        e = (ESIMMON*) malloc (sizeof(ESIMMON));
        e->rundir = strdup(".");
    }
#endif

    // At this point outpath may include eSiMon rundir ex:outpath/run1
    if (use_esimmon == 0){
        outpath = strdup(shotdir);
    }else{
        outpath = strdup(e->rundir);
    }

    // Variables to list the directory
    // common
    int simulation_ended = 0;
    char * finishfile;
    // 1D vars
    int first_check = 0;
    oneddiagfile = malloc(strlen(inpath)+strlen("/xgc.oneddiag.bp")+1);
    strcpy(oneddiagfile,inpath);
    strcat(oneddiagfile,"/xgc.oneddiag.bp");

    // 2D vars
    int next_time_step_exists = 0;
    int current_time_step = 1;
    int latest_time_step_processed = 1;
    char current_time_step_string[32];
    char latest_time_step_string[32];
    char * meshfile;
    char * datafile;
    meshfile = malloc(strlen(inpath)+strlen("/xgc.meshfile.bp")+1);
    strcpy(meshfile,inpath);
    strcat(meshfile,"/xgc.mesh.bp");
    finishfile = malloc(strlen(inpath)+strlen("/finish.sim")+1);
    strcpy(finishfile,inpath);
    strcat(finishfile,"/finish.sim");
    datafile = malloc(strlen(inpath)+strlen("/xgc.3d.00000.bp")+1);

    last_step_processed = 0;

    // while simulation ended is not true
    while (simulation_ended == 0){
        // Check if the directory exists
        if (!is_dir(inpath)){
            printf("This directory does not exist.\n");
        }
        /* Check if simulation has ended */
           if (file_exists(finishfile)){
           printf("simulation has ended %s.\n", finishfile);
           simulation_ended = 1;
        }
        // Check if meshfile exists
        if (!file_exists(meshfile)){
            printf("mesh file does not exist.\n");
            // Sleep at end of loop instead - sleep(10);
            continue;
        }

        /* 1D main body of work */
        // Check file in question 
        // last_step_processed = 7;
        if (first_check == 0){
            // Check if file exists
            if (!file_exists(oneddiagfile)){
                printf("File 'xgc.oneddiag.bp' does not exists yet in %s, waiting...\n",inpath);
                // Sleep at end of loop instead - sleep(10);
                continue;
            }else{
                // If it exists, get the original file size
                stat(oneddiagfile, &fileattrib);
                previous_file_size = fileattrib.st_size;
                first_check = 1;
            }
        }else{
            checkFile();
            last_step_processed = last_step_processed + number_steps_processed;
            number_steps_processed = 0;
            printf("last time step processed: %d\n",last_step_processed);
            printf("Waiting ...\n");
            // Sleep at end of loop instead - sleep(10);
        }

        /* 2D main body of work */
        snprintf(current_time_step_string, 32, "%0*d", 5, current_time_step);
        snprintf(latest_time_step_string, 32, "%0*d", 5, current_time_step-1);
        char * currentfile = 0; // = '\0';
        strcpy(datafile,inpath);
        strcat(datafile,"/");
        next_time_step_exists = find_next_step(current_time_step_string, &currentfile);

        // If next time step is there it is ok to open meshfile???
        // We will stay one step behind the simulation for now
        if (next_time_step_exists == 1){
            if (current_time_step == 1){
                // open mesh file
                storeMesh(meshfile);
                strcat(datafile,currentfile);
                timeFix2 = current_time_step - 1;
                plotAll2DVarsForTimeStep(datafile);
                latest_time_step_processed = timeFix2;
                current_time_step++;
            }else{
                strcat(datafile,currentfile);
                timeFix2 = current_time_step - 1;
                // Here read nphi, that will help not doing tricks to know which
                // variables to plot while waiting for schema solution 
                // Scalars do not need to be read it, we get them from the metadata 
                printf("plot all 2d vars for timestep: %s\n", datafile);
                plotAll2DVarsForTimeStep(datafile);
                latest_time_step_processed = timeFix2;
                current_time_step++;
                continue;
            }
        }else{
            //find_next_step(current_time_step_string, &currentfile);
            find_next_step(latest_time_step_string, &currentfile);
            strcat(datafile,currentfile);
            if (latest_time_step_processed < current_time_step -1){
                timeFix2 = current_time_step - 1;
                plotAll2DVarsForTimeStep(datafile);
                latest_time_step_processed++;
                printf("Waiting for new data... \n");
                // Sleep at end of loop instead - sleep(10);
            }else{
                printf("Waiting for new data... \n");
                // Sleep at end of loop instead - sleep(10);
            }
        }
        printf("Sleep here.\n");
        sleep(10);
    }

#ifdef ESIMON
   if (use_esimmon) {
       err = esimmon_finalizeRun(e);
       if (err) {
           fprintf(stderr, "ERROR at finalizeRun: errno=%d\n%s\n",
                   esimmon_errno, esimmon_get_last_errmsg());
       }

       esimmon_finalize();
   } else {
       free(e);
   }
#endif

   // List directories created 
   DIR *dp;
   struct dirent *ep;
   dp = opendir (outpath);
   char outputdir[256];
   char outputdir1[256];
   char outputdir2[256];
   int dir_has_subs = 0;
   char syscommand[4096];
   if (dp != NULL)
   {
       while (ep = readdir (dp)){
           DIR *dp2;
           struct dirent *ep2;
           if (!strcmp(ep->d_name,".") || !strcmp(ep->d_name,".."))
               continue;
           snprintf(outputdir, sizeof(outputdir), "%s/%s",outpath,ep->d_name);

           if(is_dir(outputdir)){
               dp2 = opendir(outputdir);
               dir_has_subs = 0;
               if (dp2 != NULL){
                   while (ep2 = readdir(dp2)){
                       if (!strcmp(ep2->d_name,".") || !strcmp(ep2->d_name,".."))
                           continue;
                       snprintf(outputdir1, sizeof(outputdir1), "%s/%s",outpath,ep->d_name);
                       snprintf(outputdir2, sizeof(outputdir2), "%s/%s/%s",outpath,ep->d_name,ep2->d_name);
                       if (is_dir(outputdir2)){
                           dir_has_subs = 1;
                           snprintf(syscommand, sizeof(syscommand), "/ccs/proj/e2e/rbarreto/xgc_plotter/calculator_movie.sh -d %s/%s/%s %s", outpath, ep->d_name, ep2->d_name, ep2->d_name);
                           system(syscommand);
                       }else{
                           //printf("outputdir2 is NOT a dir: %s, run it in outputdir1: %s\n",outputdir2, outputdir1);
                           break;
                       }
                   }
                   if (dir_has_subs == 0){
                       snprintf(syscommand, sizeof(syscommand), "/ccs/proj/e2e/rbarreto/xgc_plotter/calculator_movie.sh -d %s/%s %s", outpath, ep->d_name, ep->d_name);
                       system(syscommand);
                   }
                   (void) closedir(dp2);
               }
           }
       }
       (void) closedir (dp);
   }
   else
       perror ("Couldn't open the directory");

   /* Free allocated memory */
   myfree(prgname);
   myfree(outpath);
   myfree(inpath);
   myfree(xvar);
   myfree(yvar);
   return 0;
}

int find_next_step_string(char *str, char *orig)
{
  // static char buffer[4096];
  char *p;

  if(!(p = strstr(str, orig)))  // Is 'orig' even in 'str'?
    return 0;

  return 1;

}

int find_next_step( char *str, char ** returnstr){
    DIR *dp;
    struct dirent *ep;
    //dp = opendir ("../xgc_output");
    dp = opendir (inpath);
    int foundit = 0;
    if (dp != NULL)
    {
        while (ep = readdir (dp)){
            //puts (ep->d_name);
            // returnstr = strdup(ep->d_name); 
            foundit = find_next_step_string(ep->d_name, str);
            //printf("foudid is: %d\n", foundit);
            if (foundit == 1){
                *returnstr = malloc (strlen("xgc.3d.00000.bp")+1);
                strcpy(*returnstr, ep->d_name); 
                //returnstr = strdup(ep->d_name); 
                break;
            }else{
                continue;
            }
        }
        (void) closedir (dp);
    }
    else
        perror ("Couldn't open the directory");
    return foundit;
}

static int str_to_int(char *str, char *complete_str) {
    int value;
    // convert item to number
    errno = 0;
    value = strtol(str, (char **)NULL, 0);
    if (errno) {
        //printf("\n");
        fprintf(stderr, "Error: could not convert field into a value: %s from \"%s\"\n", str, complete_str);
        exit(200);
    }
    return value;
}

// parse a string "0, t3; 027" into an integer array
// of [0,3,27] and return tDimidx as index of "t" dimension in the array
// exits if parsing failes
static void parseDimSpec(char *str, int *dims, bool allow_time,
                        long *tDimidx, long *tValue, long *tStepping)
{
    char *token, *saveptr;
    char *steppingstr;
    char *s;  // copy of s; strtok modifies the string
    int  i=0;

    if (tDimidx != NULL) *tDimidx = -1;
    s = strndup(str, 256);
    token = strtok_r(s, " ,;x\t\n", &saveptr);
    while (token != NULL && i < MAX_DIMS) {
        //printf("\t|%s|", token);

        // check for time dimension
        if (token[0] == 't') {
            if (*tDimidx != -1) {
                fprintf(stderr, "Error: only one time dimension can be defined in \"%s\"\n", str);
                exit(200);
            }
            if (!allow_time) {
                fprintf(stderr, "Error: time specification is not allowed for this dimension specification: \"%s\"\n", str);
                exit(200);
            }
            *tDimidx = i;
            token++;
            //printf("TIME");
            // look for : separator for stepping info
            steppingstr = strchr(token, ':');
            if (steppingstr != NULL) {
                steppingstr[0] = '\0';
                // this puts \0 into token too! Token should contain the time value now
                if (steppingstr == NULL) {
                    fprintf(stderr, "Error: time dimension with stepping is ill defined in \"%s\". A number is needed too after the : like t0:2\n", str);
                    exit(200);
                }
                *tStepping = str_to_int(steppingstr, str);
            }
            if (token == NULL) {
                //printf("\n");
                fprintf(stderr, "Error: time dimension is ill defined in \"%s\". A number is needed too like t0\n", str);
                exit(200);
            }
            // finally get the value
            *tValue = str_to_int(token, str);
            dims[i] = *tValue;

        } else {
            // just a number
            dims[i] = str_to_int(token, str);
        }

        // get next item
        token = strtok_r(NULL, " ,;x\t\n", &saveptr);
        i++;
    }
    //if (i>0) printf("\n");

    // check if number of dims specified is larger than we can handle
    if (token != NULL) {
        fprintf(stderr, "Error: More dimensions specified in \"%s\" than we can handle (%d)\n", str, MAX_DIMS);
        exit(200);
    }
}

// parse a string "0, 255; 127" into an int[3] array
// exit if string does not conform (have exactly three values)
static void parseRGB(char *str, int rgb[3])
{
    char *token, *saveptr;
    char *s;  // copy of s; strtok modifies the string
    int  i=0;

    s = strndup(str, 256);
    token = strtok_r(s, " ,;x\t\n", &saveptr);
    //printf("Parse RGB string %s\n", s);
    while (token != NULL && i < 3) {
        //printf("\t|%s|", token);
        rgb[i] = str_to_int(token, str); // will exit if fails
        // sanity check
        if (rgb[i] != CUT_TO_BYTE(rgb[i])) {
            fprintf(stderr, "Warning: RGB value %d in \"%s\" is out of range 0..255. Will use %d\n",
                    rgb[i], str, CUT_TO_BYTE(rgb[i]));
        }
        // get next item
        token = strtok_r(NULL, " ,;x\t\n", &saveptr);
        i++;
    }
    //if (i>0) printf("\n");
    // error checks
    if (i != 3) {
        fprintf(stderr, "Error: Less than 3 values are specified in \"%s\"\n", str);
        exit(200);
    }
    if (token != NULL) {
        fprintf(stderr, "Error: More than 3 values are specified in \"%s\"\n", str);
        exit(200);
    }
}
