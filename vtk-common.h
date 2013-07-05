#ifndef __VTKCOMMON_H__
#define __VTKCOMMON_H__

extern "C" {
#   include "array.h"
}

#include "vtkDataArray.h"

#ifdef __VTKMAIN__
#   define EXT 
#else
#   define EXT extern
#endif

// Variables global for vtk-* functions
EXT int imgSizeX;
EXT int imgSizeY;
EXT double foregroundRGB[3];
EXT double backgroundRGB[3];
/*
EXT int use_polar_coordinates;    // 0: no, 1: x is radius, y is phi, other: x is phi, y is radius
EXT bool make_square_plot;        // scale actor to have equal size of x and y axis
EXT bool show_axes;               // show axes
*/

// Types
enum Colormaps { cmap_RedBlue, cmap_BlueRed, cmap_GrayScale, cmap_XGC, cmap_XGCLog, cmap_HotDesaturated };

/*************/
/* Functions */
/*************/
int getVTKType(ArrayType atype);

// Transform the Array into a vtkDataArray, types are handles, use result->GetDataType() to get it
vtkDataArray * arrayToVTKArray(Array a);

// Create and return the colormap selected by the user (or default)
//  colormap may depend on the data range, can be different for all-positive, all-negative  data
vtkLookupTable * getSelectedLookupTable( double datamin, double datamax);

// Set for Off-screen rendering. Call from each drawing class before drawing
void vtk_set_offscreen(); 

#endif
