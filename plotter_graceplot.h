#ifndef __GRACEPLOT_H__
#define __GRACEPLOT_H__

#include "common.h"

/*
void grace_plot (char *outputpath, char *title, char *subtitle, char *xname, char *yname,
                 int npoints, Array x, Array y, char* timestr,
                 double xmin, double xmax, double ymin, double ymax,
		 bool save_to_agr);
*/

/* For one image, call init() once, as many draw() as many data sets are 
 * to be plotted and save() once.
*/
 
void grace_init_plot (char *title, char *subtitle, char *xname, char *yname, float legendwidth);
void grace_draw_plot (int npoints, Array x, Array y, char *yname, char *timestr, char *legend,
                 double xmin, double xmax, double ymin, double ymax);
void grace_save_plot (char *outputpath, bool save_to_agr);

/* OPTIONAL FUNCTIONS */
// set image size
void grace_setImageSize(int x, int y);
// set foreground color with RGB (0..255 values for each color component)
void grace_setForeground( int rgb[3] );
// set background color with RGB (0..255 values for each color component)
void grace_setBackground( int rgb[3] );

#endif
