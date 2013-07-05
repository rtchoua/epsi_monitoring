#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <grace_np.h>
#ifndef __USE_ISOC99
#  define __USE_ISOC99
#endif
#include <math.h>         // ISO_C99 defines INFINITY
#include <sys/stat.h>
#include <sys/types.h>

#include "common.h"
#include "array.h"
//#include "graceplot.h"

// global vars for this module
static bool grace_started = false;
static int  grimgSizeX = 360;
static int  grimgSizeY = 300;
static int  grForegroundRGB[3] = {0, 0, 0}; 
static int  grBackgroundRGB[3] = {255, 255, 255}; 


void grace_error_function(const char *msg) {
    fprintf(stderr, "Error: xmgrace message: \"%s\"\n", msg);
}

void grace_stop() { 
    if (grace_started) {
        GraceFlush();
        GraceClose(); 
        grace_started = false;
    }
}

void grace_start() {
    //defaults write com.apple.x11 wm_ffm true
    GraceRegisterErrorFunction(grace_error_function);

    /* Start Grace with a buffer size of 20480 and open the pipe */
    if (GraceOpenVA("gracebat",20480,"-nosafe","-noask",NULL) == -1) {
        fprintf(stderr, "Can't run Grace. \n");
        exit(EXIT_FAILURE);
    }
    grace_started = true;
    atexit(grace_stop);   // ensure that we stop grace when exiting from the main program
}


void grace_setImageSize(int x, int y) {
    grimgSizeX = x;
    grimgSizeY = y;
}

void grace_setForeground( int rgb[3] ) {
    grForegroundRGB[0] = CUT_TO_BYTE(rgb[0]);
    grForegroundRGB[1] = CUT_TO_BYTE(rgb[1]);
    grForegroundRGB[2] = CUT_TO_BYTE(rgb[2]);
}

void grace_setBackground( int rgb[3] ) {
    grBackgroundRGB[0] = CUT_TO_BYTE(rgb[0]);
    grBackgroundRGB[1] = CUT_TO_BYTE(rgb[1]);
    grBackgroundRGB[2] = CUT_TO_BYTE(rgb[2]);
}
#if 0
void grace_plot (char *outputpath, char *title, char *subtitle, char *xname, char *yname,  
                 int npoints, Array x, Array y, char *timestr,
                 double xmin, double xmax, double ymin, double ymax,
                 bool save_to_agr) {
    int i;
    double yedge, xedge; // 2% of the ymax-ymin, xmax-xmin, resp.
    int xxlen=npoints;
    char xstr[128], ystr[128];

    if (!grace_started) grace_start();
    
    // If the whole data is undefined, do not produce image (would be empty)
    if ( Array_is_data_undef_all(y) ) {
        printf("All data is undef for this variable: %s %s\n", yname, timestr);
        return;
    }

    // now just plot these in grace!
    yedge = (ymax-ymin)*0.02;
    xedge = (xmax-xmin)*0.02;
    // if function is constant, yedge is zero and plot would be empty
    if (yedge == 0) yedge = ymax/2.0 + 0.55;
    //printf("\n   x: %lf - %lf\n   y: %lf - %lf\n   xedge: %lf\n   yedge: %lf\n",
    //       xmin, xmax, ymin, ymax, xedge, yedge);


    /* Send some initialization commands to Grace */


    GracePrintf("device \"PNG\" OP \"compression:9\"");
    //GracePrintf("device \"JPEG\" OP \"quality:50\"");
    GracePrintf("hardcopy device \"PNG\"");
    //GracePrintf("hardcopy device \"JPEG\"");
    //GracePrintf("page size 360,300"); 
    GracePrintf("page size %d,%d", grimgSizeX, grimgSizeY); 
    /*     GracePrintf("xaxis tick off"); */
    /*     GracePrintf("yaxis tick off"); */
    /*     GracePrintf("xaxis ticklabel off"); */
    /*     GracePrintf("yaxis ticklabel off"); */
    /*GracePrintf("view 0.137821, 0.051282, 1.157051, 0.868590");*/
    GracePrintf("view 0.15, 0.051282, 1.157051, 0.868590");
    /*     GracePrintf("s0 on"); */
    /*     GracePrintf("s0 symbol 1"); */
    /*     GracePrintf("s0 symbol size 0.3"); */
    /*     GracePrintf("s0 symbol fill pattern 1"); */
    /*     GracePrintf("s1 on"); */
    /*     GracePrintf("s1 symbol 1"); */
    /*     GracePrintf("s1 symbol size 0.3"); */
    /*     GracePrintf("s1 symbol fill pattern 1"); */
    GracePrintf("title size 1.4");
    GracePrintf("xaxis  ticklabel char size 1.36");
    GracePrintf("yaxis  ticklabel char size 1.36");
    GracePrintf("title \"%s\"",title);
    if (subtitle)
        GracePrintf("subtitle \"%s\"",subtitle);
    else
        GracePrintf("subtitle \"T=%s, [%6.2le:%6.2le],[%6.2le:%6.2le]\"",timestr,xmin,xmax,ymin,ymax);
    //GracePrintf("xaxis label \"%s\"",xname);
    //GracePrintf("yaxis label \"%s\"",yname);
    GracePrintf("map color 0 to (%d,%d,%d), \"background\"", 
                grBackgroundRGB[0], grBackgroundRGB[1], grBackgroundRGB[2]);
    GracePrintf("map color 1 to (%d,%d,%d), \"foreground\"", 
                grForegroundRGB[0], grForegroundRGB[1], grForegroundRGB[2]);
    GracePrintf("background color 0");


    for (i = 0; i < xxlen && GraceIsOpen(); i++) {
        Array_strvalue(x, i, xstr);
        Array_strvalue(y, i, ystr);
        GracePrintf("g0.s0 point %s, %s",xstr,ystr);
    }
    //GracePrintf("autoscale");
    GracePrintf("world xmin %le", xmin-xedge);
    GracePrintf("world xmax %le", xmax+xedge);
    GracePrintf("world ymin %le", ymin-yedge);
    GracePrintf("world ymax %le", ymax+yedge);
    GracePrintf("autoticks");
    GracePrintf("s0 line linewidth 4.0");
    GracePrintf("xaxis  tick major linewidth 3.0");
    GracePrintf("yaxis  tick major linewidth 3.0");
    GracePrintf("frame linewidth 2.0");
    GracePrintf("xaxis bar linewidth 3.0");
    GracePrintf("redraw");
    //if (save_to_png) {
        GracePrintf("print to \"%s.png\"",outputpath);
        GracePrintf("print");
        if (print_provenance)
            printf("--ProvenanceInfo name=\"%s.png\" size=0 date=0\n", outputpath);
    //}
    if (GraceIsOpen() && save_to_agr) {
        /* Tell Grace to save the data */
        GracePrintf("saveall \"%s.agr\"",outputpath);
    }
    GracePrintf("kill g0.s0");

}
#endif

/***************************************************

    EXPERIMENT FOR MULTI-PLOTS


****************************************************/

// variables initialized in init_plot and updated in draw_plot calls until save_plot
static bool grSubtitleGiven;
static double grXMin;
static double grXMax;
static double grYMin;
static double grYMax;
static int    grSets; // number of graphs (point sets)
static int    grInited; // true: init was called (maybe no draws afterwards)

void grace_init_plot (char *title, char *subtitle, char *xname, char *yname, float legendwidth) 
{
    float boxpos_right; // right view pos of the graph frame
    //printf("grace_init_plot: %s\n", title);
    if (!grace_started) grace_start();
    grXMin = INFINITY;
    grXMax = -INFINITY;
    grYMin = INFINITY;
    grYMax = -INFINITY;
    grSets = 0; 

    /* Send some initialization commands to Grace */
    GracePrintf("device \"PNG\" OP \"compression:9\"");
    //GracePrintf("device \"JPEG\" OP \"quality:50\"");
    GracePrintf("hardcopy device \"PNG\"");
    //GracePrintf("hardcopy device \"JPEG\"");
    //GracePrintf("page size 360,300"); 
    GracePrintf("page size %d,%d", grimgSizeX, grimgSizeY); 
    /*     GracePrintf("xaxis tick off"); */
    /*     GracePrintf("yaxis tick off"); */
    /*     GracePrintf("xaxis ticklabel off"); */
    /*     GracePrintf("yaxis ticklabel off"); */
    /*GracePrintf("view 0.137821, 0.051282, 1.157051, 0.868590");*/
    boxpos_right = 1.157051;
    if (legendwidth > 0.0) {
        GracePrintf("view 0.15, 0.051282, %f, 0.868590", boxpos_right-legendwidth);
        GracePrintf("legend on");
        GracePrintf("legend loctype view");
        GracePrintf("legend x1 %f", boxpos_right-legendwidth+0.02);
        GracePrintf("legend y1 0.86");
    } else {
        GracePrintf("view 0.15, 0.051282, %f, 0.868590", boxpos_right);
    }
    /*     GracePrintf("s0 on"); */
    /*     GracePrintf("s0 symbol 1"); */
    /*     GracePrintf("s0 symbol size 0.3"); */
    /*     GracePrintf("s0 symbol fill pattern 1"); */
    /*     GracePrintf("s1 on"); */
    /*     GracePrintf("s1 symbol 1"); */
    /*     GracePrintf("s1 symbol size 0.3"); */
    /*     GracePrintf("s1 symbol fill pattern 1"); */
    GracePrintf("title size 1.4");
    GracePrintf("xaxis  ticklabel char size 1.36");
    GracePrintf("yaxis  ticklabel char size 1.36");
    GracePrintf("title \"%s\"",title);
    if (subtitle) {
        GracePrintf("subtitle \"%s\"",subtitle);
        grSubtitleGiven = true;
    } else {
        grSubtitleGiven = false;
    }
    //GracePrintf("xaxis label \"%s\"",xname);
    //GracePrintf("yaxis label \"%s\"",yname);
    GracePrintf("map color 0 to (%d,%d,%d), \"background\"", 
                grBackgroundRGB[0], grBackgroundRGB[1], grBackgroundRGB[2]);
    GracePrintf("map color 1 to (%d,%d,%d), \"foreground\"", 
                grForegroundRGB[0], grForegroundRGB[1], grForegroundRGB[2]);
    GracePrintf("background color 0");

    GracePrintf("xaxis  tick major linewidth 3.0");
    GracePrintf("yaxis  tick major linewidth 3.0");
    GracePrintf("frame linewidth 2.0");
    GracePrintf("xaxis bar linewidth 3.0");

    grInited = true;
}


void grace_draw_plot (int npoints, Array x, Array y, char *yname, char *timestr, char *legend,
                     double xmin, double xmax, double ymin, double ymax) 
{
    int i;
    int xxlen=npoints;
    char xstr[128], ystr[128];

    // Note grSets is now the index of this plot set. 
    // Will be incremented at the end to indicate the number of sets drawn
    
    if (!grInited) {
        fprintf (stderr, "Error: grace draw was called before init\n");
        return;
    }
    // printf("grace_draw_plot %d: %s\n", grSets, yname);
    // printf("   x: %lf - %lf\n   y: %lf - %lf\n",
    //       xmin, xmax, ymin, ymax);
    
    // If the whole data is undefined, do not produce image (would be empty)
    if ( Array_is_data_undef_all(y) ) {
        printf("All data is undef for this variable: %s %s\n", yname, timestr);
        return;
    }

    // recalc min/max
    if (xmin < grXMin) grXMin = xmin;
    if (xmax > grXMax) grXMax = xmax;
    if (ymin < grYMin) grYMin = ymin;
    if (ymax > grYMax) grYMax = ymax;

    if (!grSubtitleGiven)
        GracePrintf("subtitle \"T=%s, [%6.2le:%6.2le],[%6.2le:%6.2le]\"",timestr,xmin,xmax,ymin,ymax);

    for (i = 0; i < xxlen && GraceIsOpen(); i++) {
        Array_strvalue(x, i, xstr);
        Array_strvalue(y, i, ystr);
        GracePrintf("g0.s%d point %s, %s",grSets,xstr,ystr);
    }
    GracePrintf("s%d line linewidth 4.0", grSets);
    GracePrintf("s%d line color %d", grSets, grSets%15+1);
    if (legend != NULL)
        GracePrintf("s%d legend \"%s\"", grSets, legend);

    //GracePrintf("redraw");
    //GracePrintf("print to \"i%d.png\"",grSets);
    //GracePrintf("print");

    grSets++;
}

void grace_save_plot (char *outputpath, bool save_to_agr) 
{
    int i;
    double yedge, xedge; // 2% of the ymax-ymin, xmax-xmin, resp.
    if (!grInited) 
        return;  // silently, so no warning if only 2D plots were done in main

    //printf("grace_save_plot: %s\n", outputpath);
    if (grSets > 0) {
        // draw the graphs now (if we have any data)
        // add 2% edge to the min/max values 
        yedge = (grYMax-grYMin)*0.02;
        xedge = (grXMax-grXMin)*0.02;
        // if function is constant, yedge is zero and plot would be empty
        if (yedge == 0) yedge = grYMax/2.0 + 0.55;
        //printf("\n   x: %lf - %lf\n   y: %lf - %lf\n   xedge: %lf\n   yedge: %lf\n",
        //       grXMin, grXMax, grYMin, grYMax, xedge, yedge);

        //GracePrintf("autoscale");
        GracePrintf("world xmin %le", grXMin-xedge);
        GracePrintf("world xmax %le", grXMax+xedge);
        GracePrintf("world ymin %le", grYMin-yedge);
        GracePrintf("world ymax %le", grYMax+yedge);
        GracePrintf("autoticks");
        GracePrintf("redraw");
    }
    //if (save_to_png) {
        GracePrintf("print to \"%s.png\"",outputpath);
        GracePrintf("print");
        if (print_provenance)
            printf("--ProvenanceInfo name=\"%s.png\" size=0 date=0\n", outputpath);
    //}
    if (GraceIsOpen() && save_to_agr) {
        /* Tell Grace to save the data */
        GracePrintf("saveall \"%s.agr\"",outputpath);
    }
    // clear plot (for next init, draws, save cycle)
    for (i=0; i<grSets; i++) {
        GracePrintf("kill g0.s%d", i);
    }

    grInited = false;
}




