#ifndef __VTKGRAPH_H__
#define __VTKGRAPH_H__

#ifdef __cplusplus
extern "C" {
#endif
    
#include "array.h"

typedef struct {
    char *title;               // title to put on the image
    char *xname;               // label of the X axis if appropriate
    char *yname;               // label of the Y axis if appropriate
    char *imagefilename;       // output image file 
    bool do_contour;           // do a contour plot instead of pseudo-colormap  
    int  numContours;          // number of contour levels to draw
    bool do_square;            // make the plot square by scaling x/y appropriately
    bool radius_is_first_dim;  // for polar plots, true: 1st dim of v is radius and 2nd is angle
    bool show_axes;            // show the axes: for array and rectliniear plots only 
    bool auto_zoom;            // zoom a bit to fill image. Not used if axes are plotted. 
                               //   for array, rectilinear and structured grid plots only
    bool balance_colormap;     // balance -/+ min/max values to have 0 in the middle
                               //   but only if difference is < 50% of total distance
    bool node_centered;        // handle 2D data as node centered dataset instead of
                               //  cell centered one. Makes gradient coloring 
                               //  Only for plain 2D array plotting. Others are
                               //  node centered anyway.
    bool transposed;           //  The input 2D array has been transposed before plotting
} vtk_options;
/** to initialize a vtk_options variable, instead of 
         vtk_options opts; 
    use
         vtk_options opts = VTK_OPTIONS_NEW; 
*/
#define VTK_OPTIONS_NEW { NULL, NULL, NULL, NULL, false, 5, false, true, false, true, false, true }

// init should be called before using vtk functions
void vtk_init();

// functions for settings
void vtk_setImageSize(int x, int y);
void vtk_setForeground( int rgb[3] );            // 0..255 values
void vtk_setBackground( int rgb[3] );            // 0..255 values
char * vtk_useColormap( char *colormap_name );   // redblue, bluered, gray, xgc, ret NULL on error
void vtk_usePolarCoordinates(bool x_is_phi);     // x,y 1D arrays as (radius, phi) or (phi, radius)
void vtk_makeSquarePlot(void);                   // make horizontal and vertical sizes equal
void vtk_printSettings(void);                    // for printSettings()

/*******************************/
/*  functions for 2D graphics  */
/*******************************/
/** Plot a 2D array (as structured points) 
 *  v: 2D array of values
 *  x,y: 1D arrays of X-Y Axes values (if opts.show_axes == true)
 */
int vtk_plotArray(  Array v, double min_v, double max_v, 
                    Array x, double min_x, double max_x,
                    Array y, double min_y, double max_y,
                    vtk_options opts);

/** Plot a 2D array on a rectilinear grid
 *  v: 2D array of values
 *  x: 1D array of X co-ordinates (monotonic)
 *  y: 1D array of Y co-ordinates (monotonic)
 *  1st dim of v = length of x
 *  2nd dim of v = length of y
 */  
int vtk_plotRectilinearGrid( Array v, double min_v, double max_v, 
                             Array x, double min_x, double max_x,
                             Array y, double min_y, double max_y,
                             vtk_options opts);

/** Plot a 2D array on a structured grid.
 *  v: 2D array of values
 *  x: 2D array of X co-ordinates 
 *  y: 2D array of Y co-ordinates
 *  dimensions of v = dimensions of x = dimensions of y
 */  
int vtk_plotStructuredGrid( Array v, double min_v, double max_v, 
                            Array x, double min_x, double max_x,
                            Array y, double min_y, double max_y,
                            vtk_options opts);

/** Plot a 1D array on a structured grid.
 *  v: 1D array of values
 *  x: 1D array of X co-ordinates 
 *  y: 1D array of Y co-ordinates
 *  dimensions of v = dimensions of x = dimensions of y
 */  
int vtk_plotStructuredGrid1D( Array v, double min_v, double max_v, 
                            Array x, double min_x, double max_x,
                            Array y, double min_y, double max_y,
                            vtk_options opts);

/** Plot a 2D array on a polar grid.
 *  v:   2D array of values
 *  r:   1D array of radius data
 *  phi: 1D array of angle data
 *  opt.radius_is_first_dim -> 1st dim of v = length of r and  2nd dim of v = length of phi 
 *  otherwise                  2nd dim of v = length of r and  1st dim of v = length of phi 
 */  
int vtk_plotPolarGrid( Array v,   double min_v,   double max_v,
                       Array r,   double min_r,   double max_r,
                       Array phi, double min_phi, double max_phi,
                       vtk_options opts);

/** Plot a 1D array on a triangle mesh.
 *  if length(v) is N, 
 *     coords is a Nx2 or Nx3 array of point x-y[-z] co-ordinates
 *     connections is an int32 or int16 Mx3 array of triangles (with point indexes (to coords array)) 
 */ 
int vtk_plotTriangleMesh( Array v,           double min_v,           double max_v,
                          Array coords,      double min_coords,      double max_coords,
                          Array connections, double min_connections, double max_connections,
                          vtk_options opts);


/** reset everything before making a new kind of plot.
 *  Resources are reused for the same kind of plot, so call this only if you make
 *  two different plots. 
 */
void vtk_reset();

/** debug function: check if sizeof(char) in C equals sizeof(bool) in C++. */
void error_check_sizeofbool( int size_of_char_in_C );

#ifdef __cplusplus
}
#endif

#endif
