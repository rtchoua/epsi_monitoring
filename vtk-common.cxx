#define __VTKMAIN__
#include "vtk-common.h"
#include "vtk-graph.h"

extern "C" {
#   include "common.h" // verbose 
#   include "array.h"
}

#include "vtkType.h"
#include "vtkDataArray.h"
#include "vtkTypeInt16Array.h"
#include "vtkTypeInt32Array.h"
#include "vtkTypeInt64Array.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkGraphicsFactory.h"
#include "vtkImagingFactory.h"
#include "vtkLookupTable.h"


#include <stdio.h>
#include <stdint.h>
#include <regex.h>


void vtk_init() {
    // set default values
    imgSizeX = 518;
    imgSizeY = 518;
    foregroundRGB[0] = 0.9;
    foregroundRGB[1] = 0.9;
    foregroundRGB[2] = 0.9;
    backgroundRGB[0] = 0.9;
    backgroundRGB[1] = 0.9;
    backgroundRGB[2] = 0.9;
    /*
    numContours = 5;
    use_polar_coordinates = 0;
    make_square_plot = false;
    show_axes = false;
    */
}

void vtk_set_offscreen() {
    static bool calledAlready = false; 

    if (!calledAlready) {
        // Graphics Factory
        vtkGraphicsFactory * graphics_factory
                = vtkGraphicsFactory::New();
        graphics_factory->SetOffScreenOnlyMode( 1);
        graphics_factory->SetUseMesaClasses( 1 );
    
        // Imaging Factory
        vtkImagingFactory * imaging_factory 
                = vtkImagingFactory::New();
        imaging_factory->SetUseMesaClasses( 1 );

        calledAlready = true;
    }
}

void vtk_setImageSize(int x, int y) {
    imgSizeX = x;
    imgSizeY = y;
}

void vtk_setForeground( int rgb[3] ) {
    foregroundRGB[0] = CUT_TO_BYTE(rgb[0])/255.0; // 0.0 - 1.0 value for VTK
    foregroundRGB[1] = CUT_TO_BYTE(rgb[1])/255.0;
    foregroundRGB[2] = CUT_TO_BYTE(rgb[2])/255.0;
}

void vtk_setBackground( int rgb[3] ) {
    backgroundRGB[0] = CUT_TO_BYTE(rgb[0])/255.0; // 0.0 - 1.0 value for VTK
    backgroundRGB[1] = CUT_TO_BYTE(rgb[1])/255.0;
    backgroundRGB[2] = CUT_TO_BYTE(rgb[2])/255.0;
}


/*
void vtk_usePolarCoordinates(bool x_is_phi) {
    use_polar_coordinates = 1+x_is_phi; // 1: x is radius, 2: x is phi
    //printf("use_polar_coordinates = %d\n", use_polar_coordinates);
}

void vtk_makeSquarePlot() {
    make_square_plot = true;
}

void vtk_setNumberOfContours( int nContours ) {
    if (nContours < 0) {
        fprintf(stderr, "Error: number of contours should be between 0 and 64: set to 0 now\n");
        nContours = 0;
    }
    if (nContours > 64) {
        fprintf(stderr, "Error: number of contours should be between 0 and 64: set to 64 now\n");
        nContours = 64;
    }
    numContours = nContours;
}
*/

int getVTKType(ArrayType atype) {
    int vtkType;
    switch (atype) {
        case int16Array:
            vtkType = VTK_TYPE_INT16;
            break;
        case int32Array:
            vtkType = VTK_TYPE_INT32;
            break;
        case int64Array:
            vtkType = VTK_TYPE_INT64;
            break;
        case floatArray:
            vtkType = VTK_FLOAT;
            break;
        case doubleArray:
            vtkType = VTK_DOUBLE;
            break;
        default:
            fprintf(stderr, "getVTKType: Error: unknown array type (type idx=%d) \n", atype);
            exit(102);
    }
    return vtkType;
}

vtkDataArray * arrayToVTKArray(Array a) {
    // get info and data array pointer from Array
    int ndims = Array_getDimensions(a);
    int size = 1;
    for (int i=0; i<ndims; i++) size *= Array_getDimensionSize(a, i);
    ArrayType atype = Array_getType(a);
    void *data = Array_getDataPointer(a);

    // create an appropriate VTK array and set its array to the pointer
    vtkDataArray * vtkarray;
    vtkTypeInt16Array * ai16;
    vtkTypeInt32Array * ai32;
    vtkTypeInt64Array * ai64;
    vtkFloatArray * afloat;
    vtkDoubleArray * adouble;
    switch (atype) {
        case int16Array:
            ai16 = vtkTypeInt16Array::New();
            ai16->SetArray( (int16_t *)data, size, 1);
            vtkarray = ai16;
            break;
        case int32Array:
            ai32 = vtkTypeInt32Array::New();
            ai32->SetArray( (int32_t *)data, size, 1);
            vtkarray = ai32;
            break;
        case int64Array:
            ai64 = vtkTypeInt64Array::New();
            ai64->SetArray( (long long int *)data, size, 1);
            vtkarray = ai64;
            break;
        case floatArray:
            afloat = vtkFloatArray::New();
            afloat->SetArray( (float *)data, size, 1);
            vtkarray = afloat;
            break;
        case doubleArray:
            adouble = vtkDoubleArray::New();
            adouble->SetArray( (double *)data, size, 1);
            vtkarray = adouble;
            break;
        default:
            fprintf(stderr, "arrayToVTKArray: Error: unknown array type (type idx=%d) \n", atype);
            exit(202);
    }
    return vtkarray;
}


/*************************/
/*       COLORMAPS       */
/*************************/
static enum Colormaps selected_colormap = cmap_BlueRed;
struct cmapname {
    enum Colormaps cmap;
    char *name;
};

static struct cmapname cmaps[] = {
     { cmap_RedBlue,        "RedBlue"}
    ,{ cmap_BlueRed,        "BlueRed"}
    ,{ cmap_GrayScale,      "Gray"} 
    ,{ cmap_XGC,            "XGC" }
    ,{ cmap_XGCLog,         "XGCLog" }
    ,{ cmap_HotDesaturated, "HotDesaturated" }
    ,{ cmap_RedBlue,        NULL}  // The last item to stop loops on this array
};

char * vtk_useColormap( char * colormap_name ) {
    int ex, i;
    char buf[256];
    regex_t reg;
    ex = regcomp( &reg, colormap_name, REG_ICASE | REG_NOSUB);
    if (ex) {
        regerror(ex, &reg, buf, 256);
        fprintf(stderr, "Error: colormap name '%s' is invalid regular expression: %s\n", colormap_name, buf);
        return NULL;
    }
    i = 0;
    ex = 1; // not found
    while (cmaps[i].name != NULL) {
        if ( !regexec( &reg, cmaps[i].name, 0, NULL, 0) ) {
            selected_colormap = cmaps[i].cmap;
            ex = 0; // found
            break;
        }
        i++;
    }
    if (ex) { // if not found
        fprintf(stderr, "No colormap name matched the argument %s\n", colormap_name);
        selected_colormap = cmap_BlueRed;
        fprintf(stderr, "Use default %s colormap\n", cmaps[selected_colormap].name);
    }
    return cmaps[selected_colormap].name;
}

static vtkLookupTable * lut_BlueRed = NULL;
static vtkLookupTable * lut_RedBlue = NULL;
static vtkLookupTable * lut_GrayScale = NULL;
static vtkLookupTable * lut_HotDesaturated = NULL;  // VisIt hot_desaturated colormap
static vtkLookupTable * lut_XGC = NULL;
static vtkLookupTable * lut_XGC_Positive = NULL;
static vtkLookupTable * lut_XGC_Negative = NULL;
static vtkLookupTable * lut_XGCLog = NULL;
static vtkLookupTable * lut_XGCLog_Positive = NULL;
static vtkLookupTable * lut_XGCLog_Negative = NULL;

/* XGC and XGCLog colormap */
static const int XGC_colors = 9;        // number of anchor colors in XGC_RGBA
static const int XGC_numColors =64;   // how many colors to generate
static const int XGC_zero_coloridx = 4; // the middle (0.0) color XGC_RGBA[4]

static double XGC_RGBA[XGC_colors][4] = {
    {         0.0,         0.0,         0.0, 1.0 } // black
   ,{ 142.0/255.0,  41.0/255.0, 178.0/255.0, 1.0 } // magenta
   ,{   0.0/255.0,   0.0/255.0,         1.0, 1.0 } // blue
   ,{   0.0/255.0,         1.0,         1.0, 1.0 } // cyan
   ,{         1.0,         1.0,         1.0, 1.0 } // white
   ,{         1.0,         0.0,         0.0, 1.0 } // red
   ,{         1.0, 119.0/255.0,         0.0, 1.0 } // orange
   ,{         1.0,         1.0,         0.0, 1.0 } // yellow
   ,{ 167.0/255.0,         1.0,         0.0, 1.0 } // light green
   //,{         0.6,         0.9,         0.9, 1.0 } // a gray close to white (only the end target for logarithmic map)
};

/* Hot Desaturated colormap */
static const int HotDesaturated_colors = 9;        // number of anchor colors in HotDesaturated_RGBA
static const int HotDesaturated_numColors = 64;   // how many colors to generate
//static const int HotDesaturated_zero_coloridx = 4; // the middle (0.0) color HotDesaturated_RGBA[4]

static double HotDesaturated_RGBA[HotDesaturated_colors][4] = {
    {  71.0/255.0,  71.0/255.0, 219.0/255.0, 1.0 } // purplish blue
   ,{         0.0,         0.0,  91.0/255.0, 1.0 } // dark blue
   ,{         0.0,         1.0,         1.0, 1.0 } // cyan
   ,{         0.0, 127.0/255.0,         0.0, 1.0 } // dark green
   ,{         1.0,         1.0,         0.0, 1.0 } // yellow
   ,{         1.0,  96.0/255.0,         0.0, 1.0 } // dark orange
   ,{ 107.0/255.0,         0.0,         0.0, 1.0 } // dark red
   ,{ 224.0/255.0,  76.0/255.0,  76.0/255.0, 1.0 } // pinkish 
   ,{         1.0,         1.0,         1.0, 1.0 } // white
};


/** Generate a linear colormap from 'nColors' colors defined in rgbacolors array.
 *  rgbacolors should be a double[nColors][4] matrix, 
 *  rgbacolors[i] defines the RGBA values of the ith color. 
 */
static vtkLookupTable * generateLinearColormap( int nColors, double (*rgbacolors)[4], int colormapsize) 
{
    vtkLookupTable * lut;
    double rstep, gstep, bstep, astep;
    int nSectionColor = colormapsize / (nColors-1);
    //printf("linear colormap of %d colors, size %d, section size=%d\n", nColors, colormapsize, nSectionColor);
    lut = vtkLookupTable::New();
    lut->SetNumberOfColors(colormapsize);
    lut->Build();
    for (int c = 0; c < nColors-1; c++) {
        // build color transition from color c to c+1
        //printf("  color %d\n", c);
        //printf("rgba %g %g %g %g\n", rgbacolors[c][0], rgbacolors[c][1], rgbacolors[c][2],rgbacolors[c][3]);
        rstep = (rgbacolors[c+1][0] - rgbacolors[c][0]) / nSectionColor;
        gstep = (rgbacolors[c+1][1] - rgbacolors[c][1]) / nSectionColor;
        bstep = (rgbacolors[c+1][2] - rgbacolors[c][2]) / nSectionColor;
        astep = (rgbacolors[c+1][3] - rgbacolors[c][3]) / nSectionColor;
        for (int i = 0; i < nSectionColor; i++) {
            lut->SetTableValue( c*nSectionColor + i, 
                                rgbacolors[c][0] + i*rstep,
                                rgbacolors[c][1] + i*gstep,
                                rgbacolors[c][2] + i*bstep,
                                rgbacolors[c][3] + i*astep
                              );

        }
    }
    return lut;
}

/** Generate a logarithmic colormap from 'nColors' colors defined in rgbacolors array.
 *  rgbacolors should be a double[nColors][4] matrix, 
 *  rgbacolors[i] defines the RGBA values of the ith color. 
 *  centercolor is the ith color, which is the center of the logarithmic scale.
 */
static vtkLookupTable * generateLogarithmicColormap( int nColors, double (*rgbacolors)[4], int centercolor, int colormapsize) 
{
    vtkLookupTable * lut;
    double rstep, gstep, bstep, astep;
    int ratios[nColors-1];  // ratio for number of colors between two defined colors
    int ratiosum;
    int basecolornum;
    int idx;
    int nSectionColor;
    // colors [ ...   c-3                c-2        c-1    c    c+1        c+2                c+3  ...]
    // ratios [ ...            4x              2x,      1x,   1x,      2x,           4x            ...]
    if (centercolor > 0) {
        ratios[centercolor-1] = 1;
        ratiosum = 1;
        for (int i=centercolor-2; i>=0; i--) {
            ratios[i] = 2*ratios[i+1];
            ratiosum += ratios[i];
        }
    }
    ratios[centercolor] = 1;
    ratiosum++;
    for (int i=centercolor+1; i<nColors-1; i++) {
        ratios[i] = 2*ratios[i-1];
        ratiosum += ratios[i];
    }
    
    // make the number of total colors in the map a proper multiply of ratiosum
    //   total number of colors = basecolornum * ratiosum
    basecolornum = colormapsize / ratiosum; // num of colors for center value
    //printf("Ratio sum=%d  base colornum=%d  total=%d\n", ratiosum, basecolornum, basecolornum*ratiosum);
   
    lut = vtkLookupTable::New();
    lut->SetNumberOfColors(basecolornum*ratiosum);
    lut->Build();
    idx = 0;
    for (int c = 0; c < nColors-1; c++) {
        nSectionColor = basecolornum * ratios[c];
        //printf("Color=%d, section colornum=%d\n", c, nSectionColor);
        // build color transition from color c to c+1
        rstep = (rgbacolors[c+1][0] - rgbacolors[c][0]) / nSectionColor;
        gstep = (rgbacolors[c+1][1] - rgbacolors[c][1]) / nSectionColor;
        bstep = (rgbacolors[c+1][2] - rgbacolors[c][2]) / nSectionColor;
        astep = (rgbacolors[c+1][3] - rgbacolors[c][3]) / nSectionColor;
        for (int i = 0; i < nSectionColor; i++) {
            lut->SetTableValue( idx,
                                rgbacolors[c][0] + i*rstep,
                                rgbacolors[c][1] + i*gstep,
                                rgbacolors[c][2] + i*bstep,
                                rgbacolors[c][3] + i*astep
                              );
            idx++;

        }
    }
    return lut;
}

vtkLookupTable * getSelectedLookupTable( double datamin, double datamax) {
    vtkLookupTable * lut;
    double rstep, gstep, bstep, astep;
    int nSectionColor;

    switch (selected_colormap) {
        case cmap_BlueRed:
            if (lut_BlueRed == NULL) {
                lut_BlueRed = vtkLookupTable::New();
                lut_BlueRed->SetHueRange( 0.6667, 0.0 );
                lut_BlueRed->Build();
            }
            lut = lut_BlueRed;
            break;

        case cmap_RedBlue:
            if (lut_RedBlue == NULL) {
                lut_RedBlue = vtkLookupTable::New();
                lut_RedBlue->SetHueRange( 0.0, 0.6667 );  // unnecessary as this is the default
                lut_RedBlue->Build();
            }
            lut = lut_RedBlue;
            break;

        case cmap_GrayScale:
            if (lut_GrayScale == NULL) {
                lut_GrayScale = vtkLookupTable::New();
                lut_GrayScale->SetHueRange( 0.0, 0.0 );
                lut_GrayScale->SetSaturationRange( 0.0, 0.0 );
                lut_GrayScale->SetValueRange( 0.0, 1.0 );
                lut_GrayScale->Build();
            }
            lut = lut_GrayScale;
            break;

        case cmap_HotDesaturated:
            if (lut_HotDesaturated == NULL) {
                lut_HotDesaturated = generateLinearColormap( HotDesaturated_colors, HotDesaturated_RGBA, HotDesaturated_numColors);
            }
            lut = lut_HotDesaturated;
            break;

        case cmap_XGC:
            // should make different colormaps for all-positive, all-negative and pos-neg data
            
            if (datamin >= 0) { 
                // all positive: = BlueRed colormap
                
                if (lut_XGC_Positive == NULL) {
                    lut_XGC_Positive = vtkLookupTable::New();
                    lut_XGC_Positive->SetHueRange( 0.6667, 0.0 );
                    lut_XGC_Positive->Build();
                    /*
                    nSectionColor = XGC_numColors / (XGC_colors-1);
                    lut_XGC_Positive = vtkLookupTable::New();
                    lut_XGC_Positive->SetNumberOfColors((XGC_colors-XGC_zero_coloridx-1) * nSectionColor);
                    lut_XGC_Positive->Build();
                    for (int c = XGC_zero_coloridx; c < XGC_colors-1; c++) {
                        // build color transition from color c to c+1
                        rstep = (XGC_RGBA[c+1][0] - XGC_RGBA[c][0]) / nSectionColor;
                        gstep = (XGC_RGBA[c+1][1] - XGC_RGBA[c][1]) / nSectionColor;
                        bstep = (XGC_RGBA[c+1][2] - XGC_RGBA[c][2]) / nSectionColor;
                        astep = (XGC_RGBA[c+1][3] - XGC_RGBA[c][3]) / nSectionColor;
                        for (int i = 0; i < nSectionColor; i++) {
                            lut_XGC_Positive->SetTableValue( (c-XGC_zero_coloridx)*nSectionColor + i, 
                                                    XGC_RGBA[c][0] + i*rstep,
                                                    XGC_RGBA[c][1] + i*gstep,
                                                    XGC_RGBA[c][2] + i*bstep,
                                                    XGC_RGBA[c][3] + i*astep
                                                  );

                        }
                    }
                    */
                }
                lut_XGC_Positive = lut_BlueRed;
                lut = lut_XGC_Positive;
            } else if (datamax <= 0) { 
                // all negative
                if (lut_XGC_Negative == NULL) {
                    nSectionColor = XGC_numColors / (XGC_colors-1);
                    lut_XGC_Negative = vtkLookupTable::New();
                    lut_XGC_Negative->SetNumberOfColors(XGC_zero_coloridx * nSectionColor);
                    lut_XGC_Negative->Build();
                    for (int c = 0; c < XGC_zero_coloridx; c++) {
                        // build color transition from color c to c+1
                        rstep = (XGC_RGBA[c+1][0] - XGC_RGBA[c][0]) / nSectionColor;
                        gstep = (XGC_RGBA[c+1][1] - XGC_RGBA[c][1]) / nSectionColor;
                        bstep = (XGC_RGBA[c+1][2] - XGC_RGBA[c][2]) / nSectionColor;
                        astep = (XGC_RGBA[c+1][3] - XGC_RGBA[c][3]) / nSectionColor;
                        for (int i = 0; i < nSectionColor; i++) {
                            lut_XGC_Negative->SetTableValue( c*nSectionColor + i, 
                                                    XGC_RGBA[c][0] + i*rstep,
                                                    XGC_RGBA[c][1] + i*gstep,
                                                    XGC_RGBA[c][2] + i*bstep,
                                                    XGC_RGBA[c][3] + i*astep
                                                  );

                        }
                    }
                    /*
                    lut_XGC_Negative->SetHueRange( 0.6667, 0.6667 );
                    lut_XGC_Negative->SetSaturationRange( 1.0, 0.0 );
                    lut_XGC_Negative->Build();
                    */
                }
                lut = lut_XGC_Negative;
            } else  {
                // negative - positive (LINEAR SCALE COLORMAP)
                if (lut_XGC == NULL) {
                    
                    lut_XGC = generateLinearColormap( XGC_colors, XGC_RGBA, XGC_numColors);
                    /*
                    nSectionColor = XGC_numColors / (XGC_colors-1);
                    lut_XGC = vtkLookupTable::New();
                    lut_XGC->SetNumberOfColors(XGC_numColors);
                    lut_XGC->Build();
                    for (int c = 0; c < XGC_colors-1; c++) {
                        // build color transition from color c to c+1
                        rstep = (XGC_RGBA[c+1][0] - XGC_RGBA[c][0]) / nSectionColor;
                        gstep = (XGC_RGBA[c+1][1] - XGC_RGBA[c][1]) / nSectionColor;
                        bstep = (XGC_RGBA[c+1][2] - XGC_RGBA[c][2]) / nSectionColor;
                        astep = (XGC_RGBA[c+1][3] - XGC_RGBA[c][3]) / nSectionColor;
                        for (int i = 0; i < nSectionColor; i++) {
                            lut_XGC->SetTableValue( c*nSectionColor + i, 
                                                    XGC_RGBA[c][0] + i*rstep,
                                                    XGC_RGBA[c][1] + i*gstep,
                                                    XGC_RGBA[c][2] + i*bstep,
                                                    XGC_RGBA[c][3] + i*astep
                                                  );

                        }
                    }
                    */
                }
                lut = lut_XGC;
            }
            break;

        case cmap_XGCLog: /* LOGARITHMIC SCALE COLORMAP for XGC*/
            // should make different colormaps for all-positive, all-negative and pos-neg data
            
            if (datamin >= 0) { 
                // all positive, STILL LINEAR (SHOULD BE LOGARITHMIC ?)
                if (lut_XGCLog_Positive == NULL) {
                    lut_XGCLog_Positive = vtkLookupTable::New();
                    lut_XGCLog_Positive->SetHueRange( 0.6667, 0.0 );
                    lut_XGCLog_Positive->Build();
                }
                lut_XGCLog_Positive = lut_BlueRed;
                lut = lut_XGCLog_Positive;
            } else if (datamax <= 0) { 
                // all negative STILL LINEAR
                if (lut_XGCLog_Negative == NULL) {
                    nSectionColor = XGC_numColors / (XGC_colors-1);
                    lut_XGCLog_Negative = vtkLookupTable::New();
                    lut_XGCLog_Negative->SetNumberOfColors(XGC_zero_coloridx * nSectionColor);
                    lut_XGCLog_Negative->Build();
                    for (int c = 0; c < XGC_zero_coloridx; c++) {
                        // build color transition from color c to c+1
                        rstep = (XGC_RGBA[c+1][0] - XGC_RGBA[c][0]) / nSectionColor;
                        gstep = (XGC_RGBA[c+1][1] - XGC_RGBA[c][1]) / nSectionColor;
                        bstep = (XGC_RGBA[c+1][2] - XGC_RGBA[c][2]) / nSectionColor;
                        astep = (XGC_RGBA[c+1][3] - XGC_RGBA[c][3]) / nSectionColor;
                        for (int i = 0; i < nSectionColor; i++) {
                            lut_XGCLog_Negative->SetTableValue( c*nSectionColor + i, 
                                                    XGC_RGBA[c][0] + i*rstep,
                                                    XGC_RGBA[c][1] + i*gstep,
                                                    XGC_RGBA[c][2] + i*bstep,
                                                    XGC_RGBA[c][3] + i*astep
                                                  );

                        }
                    }
                }
                lut = lut_XGCLog_Negative;
            } else  {
                // negative - positive (LOGARITHMIC SCALE COLORMAP)
                if (lut_XGCLog == NULL) {
                    lut_XGCLog = generateLogarithmicColormap(XGC_colors, XGC_RGBA, XGC_zero_coloridx, XGC_numColors);
                    /*
                    int ratios[XGC_colors];
                    int ratiosum = 1;
                    int basecolornum;
                    int idx;
                    //  [ ..., 16, 8, 4, 2, 1, 2, 4, 8, 16, ...]
                    ratios[XGC_zero_coloridx] = 1;
                    for (int i=XGC_zero_coloridx-1; i>=0; i--) {
                        ratios[i] = 2*ratios[i+1];
                        ratiosum += ratios[i];
                    }
                    for (int i=XGC_zero_coloridx+1; i<XGC_colors; i++) {
                        ratios[i] = 2*ratios[i-1];
                        ratiosum += ratios[i];
                    }
                    // make the number of total colors in the map a proper multiply of ratiosum
                    //   total number of colors = basecolornum * ratiosum
                    basecolornum = XGC_numColors / ratiosum; // num of colors for center value
                    printf("Ratio sum=%d  base colornum=%d  total=%d\n", ratiosum, basecolornum, basecolornum*ratiosum);
                   
                    lut_XGCLog = vtkLookupTable::New();
                    lut_XGCLog->SetNumberOfColors(basecolornum*ratiosum);
                    lut_XGCLog->Build();
                    idx = 0;
                    for (int c = 0; c < XGC_colors; c++) {
                        nSectionColor = basecolornum * ratios[c];
                        printf("Color=%d, section colornum=%d\n", c, nSectionColor);
                        // build color transition from color c to c+1
                        rstep = (XGC_RGBA[c+1][0] - XGC_RGBA[c][0]) / nSectionColor;
                        gstep = (XGC_RGBA[c+1][1] - XGC_RGBA[c][1]) / nSectionColor;
                        bstep = (XGC_RGBA[c+1][2] - XGC_RGBA[c][2]) / nSectionColor;
                        astep = (XGC_RGBA[c+1][3] - XGC_RGBA[c][3]) / nSectionColor;
                        for (int i = 0; i < nSectionColor; i++) {
                            lut_XGCLog->SetTableValue( idx,
                                                    XGC_RGBA[c][0] + i*rstep,
                                                    XGC_RGBA[c][1] + i*gstep,
                                                    XGC_RGBA[c][2] + i*bstep,
                                                    XGC_RGBA[c][3] + i*astep
                                                  );
                            idx++;

                        }
                    }
                    printf("idx=%d\n",idx);
                    */
                }
                lut = lut_XGCLog;
            }
            break;

        default:
            fprintf(stderr, "Invalid colormap is selected\n");
            exit(201);
    }

    return lut;
}


void vtk_printSettings() {
    printf("  Colormap %s\n", cmaps[selected_colormap].name); 
}

/** debug function: check if sizeof(char) in C equals sizeof(bool) in C++. */
void error_check_sizeofbool( int size_of_char_in_C ) {
    int sbool = sizeof(bool);
    if (sbool != size_of_char_in_C) {
        fprintf(stderr, "ERROR: Size of char in C != size of bool in C++. Use compilers where this is true.\n");
        exit(-1);
    }
}
