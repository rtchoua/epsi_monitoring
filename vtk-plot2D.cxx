#include <string.h>
#include <stdint.h>

#include "vtk-graph.h"
#include "vtk-common.h"
extern "C" {
#    include "array.h"
#    include "timer.h"
}

#include "vtkActor.h"
#include "vtkActor2D.h"
#include "vtkAxisActor2D.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkContourFilter.h"
#include "vtkCubeAxesActor2D.h"
#include "vtkDataSetMapper.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkLookupTable.h"
#include "vtkMarchingSquares.h"
#include "vtkPNGWriter.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPointData.h"
#include "vtkProperty2D.h"
#include "vtkRectilinearGrid.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkScalarBarActor.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridGeometryFilter.h"
#include "vtkTextProperty.h" 
#include "vtkWindowToImageFilter.h"

static vtkDataSetMapper         * dataMapper    = NULL;
static vtkPolyDataMapper        * polyMapper    = NULL;
static vtkMarchingSquares       * contourMS     = NULL;
static vtkContourFilter         * contour       = NULL;
static vtkActor                 * dataActor     = NULL;
static vtkScalarBarActor        * scalarBar     = NULL;
static vtkRenderer              * renderer      = NULL;
static vtkRenderer              * barrenderer   = NULL;
static vtkRenderWindow          * render_window = NULL;
static vtkWindowToImageFilter   * win_2_image   = NULL;
static vtkPNGWriter             * pngwriter     = NULL;

#define MYFREE(var)   if (var!=NULL) {var->Delete(); var=NULL;}
void vtk_reset() 
{
    MYFREE (dataMapper)
    MYFREE (polyMapper)
    MYFREE (contourMS)
    MYFREE (contour)
    MYFREE (dataActor)
    MYFREE (scalarBar)
    MYFREE (renderer)
    MYFREE (barrenderer)
    MYFREE (render_window)
    MYFREE (win_2_image)
    MYFREE (pngwriter)
}


// declarations
static int vtk_plot2D( vtkMapper * mapper, vtkCubeAxesActor2D * axesActor, vtk_options opts, 
                       int dim_x,    int dim_y, 
                       double min_v, double max_v,
                       double min_x, double max_x,
                       double min_y, double max_y);

/**
 *   Plot a simple array (STRUCTURED POINTS)
 */
int vtk_plotArray(  Array v, double min_v, double max_v, 
                    Array x, double min_x, double max_x, // X-axis values only
                    Array y, double min_y, double max_y, // Y-axis values only
                    vtk_options opts)
{
    int retval, dims[2];
    for (int i=0; i<2; i++) dims[i] = Array_getDimensionSize(v, i);

    // make a vtkDataArray from the Array 
    vtkDataArray * scalars = arrayToVTKArray(v);

    // Structured Points 
    // C array dim0xdim1 --> Fortran style arrays in VTK dim1xdim0
    vtkStructuredPoints * spoints 
            = vtkStructuredPoints::New();
    if (opts.node_centered)
        spoints->SetDimensions(dims[1],dims[0],1);
    else 
        spoints->SetDimensions(dims[1]+1,dims[0]+1,1);
    //spoints->SetOrigin(-1,-1,1);
    //spoints->SetSpacing( 2.0/dims[0], 2.0/dims[1], 1);
    spoints->SetScalarType( scalars->GetDataType() );

    vtkPointData * pointData = spoints->GetPointData();
    vtkCellData * cellData = spoints->GetCellData();
    if (opts.node_centered) {
        // set point data of structured points as the vtkFloatArray
        pointData->SetScalars(scalars);
        // get range of scalars in points (use min_v and max_v now)
        //double * drange = pointData->GetScalars()->GetRange();
    } else {
        // set point data of structured points as the vtkFloatArray
        cellData->SetScalars(scalars);
    }

    vtkMapper * mapper;

    if (opts.do_contour) {
        if (contourMS == NULL) 
            contourMS = vtkMarchingSquares::New();
        contourMS->SetInput( spoints );
        contourMS->GenerateValues( opts.numContours, min_v, max_v);
        if (polyMapper == NULL)
            polyMapper = vtkPolyDataMapper::New();
        polyMapper->SetInput( contourMS->GetOutput() );
        mapper = polyMapper;
    } else {
        // Data mapper
        if (dataMapper == NULL) {
            dataMapper = vtkDataSetMapper::New();
        } 
        dataMapper->SetInput(spoints);
        mapper = dataMapper;
    }

    // create axes
    vtkCubeAxesActor2D  * axesActor = NULL;
    if (opts.show_axes) {
        // Note that axes are flipped for the arrays (2D C array -> VTK 2D array)
        double bounds[6]={0.0, dims[1]-1.0, 0, dims[0]-1.0, 0.0, 0.0};
        double ranges[6]={min_y, max_y, min_x, max_x, 0.0, 0.0};
        printf("### Axis bounds: x=%f..%f  y=%f..%f\n", bounds[0], bounds[1], bounds[2], bounds[3]);
        printf("### Axis ranges: x=%f..%f  y=%f..%f\n", min_y, max_y, min_x, max_x);
        axesActor = vtkCubeAxesActor2D::New();
        axesActor->SetBounds(bounds);
        axesActor->SetRanges(ranges);
        axesActor->SetUseRanges(1);
        axesActor->SetZAxisVisibility(0);
        if (x != NULL && y != NULL && 
            (  (Array_getType(x) != int16Array && Array_getType(x) != int32Array && Array_getType(x) != int64Array) 
            || (Array_getType(y) != int16Array && Array_getType(y) != int32Array && Array_getType(y) != int64Array)
            )
           ) {
            axesActor->SetLabelFormat("%g"); // non-integer format
        } else {
            axesActor->SetLabelFormat("%.0f"); // both axes are integers
        }

        // common settings for all axes
        vtkTextProperty * txtp = vtkTextProperty::New();
        txtp->SetColor(foregroundRGB);
        printf("Orientation = %f\n", txtp->GetOrientation());
        txtp->SetOrientation(30.0);
        axesActor->SetAxisLabelTextProperty(txtp);
        axesActor->SetAxisTitleTextProperty(txtp);
        //axesActor->SetFontFactor(1.4); // txtp->SetFontSize has no effect 
        axesActor->SetNumberOfLabels(4);
        //printf("label format = %s\n", axesActor->GetLabelFormat());

        // X-axis settings
        vtkAxisActor2D  * xAxisActor = axesActor->GetXAxisActor2D();
        xAxisActor->SetLabelFactor(1.0);
        //xAxisActor->SetLabelFormat("x%0.f");  // has no effect on the CubeAxisActor2D

        // Y-axis settings
        vtkAxisActor2D  * yAxisActor = axesActor->GetYAxisActor2D();
        yAxisActor->SetLabelFactor(1.0);
        txtp = vtkTextProperty::New();
        txtp->SetOrientation(90.0);

        // set axis names
        // Note that axis names are flipped too but only if not opts.transposed
        char *hname; // horizontal axis name
        char xfixname[] = "x";
        char yfixname[] = "y";
        if (opts.xname != NULL) {
            if (opts.transposed) hname = opts.xname;
            else                 hname = opts.yname;
        } else {
            if (opts.transposed) hname = xfixname;
            else                 hname = yfixname;
        }
        axesActor->SetXLabel(hname);

        char *vname; // vertical axis name
        if (opts.yname != NULL) {
            if (opts.transposed) vname = opts.yname;
            else                 vname = opts.xname;
        } else {
            if (opts.transposed) vname = yfixname;
            else                 vname = xfixname;
        }
        axesActor->SetYLabel(vname);
    }

    // debug
    if (verbose>1) { 
    printf("    vtk_plotArray:\n");
    printf("      v dims:           %dx%d = %d points\n", dims[0], dims[1], dims[0]*dims[1]);
    printf("      scalars:\n");
    if (opts.node_centered) 
        printf("      pointdata:        tuples = %d\n", pointData->GetAbstractArray(0)->GetNumberOfTuples());
    else
        printf("      celldata:         tuples = %d\n", cellData->GetAbstractArray(0)->GetNumberOfTuples());
    printf("      spoints: points = %d\n", spoints->GetNumberOfPoints());
    }

    scalars->Delete();
    spoints->Delete();

    retval =  vtk_plot2D( mapper, axesActor, opts, dims[0], dims[1], min_v, max_v, 
                       0.0, (double)(dims[0])-1.0, 0.0, (double)(dims[1])-1.0);

    MYFREE (axesActor);
    return retval;
}

/**
 *   Plot a 2D array over a RECTILINEAR GRID
 *   1D x and y arrays should be supplied for a 2D array.
 */
int vtk_plotRectilinearGrid( Array v, double min_v, double max_v, 
                             Array x, double min_x, double max_x,
                             Array y, double min_y, double max_y,
                             vtk_options opts)
{
    int dims[2];
    for (int i=0; i<2; i++) dims[i] = Array_getDimensionSize(v, i);

    // make a vtkDataArray from the Array 
    vtkDataArray * scalars = arrayToVTKArray(v);
    vtkDataArray * da_x = arrayToVTKArray(x);
    vtkDataArray * da_y = arrayToVTKArray(y);
    
    // Rectilinear Grid
    vtkRectilinearGrid * rgrid 
            = vtkRectilinearGrid::New();
    rgrid->SetDimensions(dims[1],dims[0],1);
    rgrid->SetXCoordinates(da_y);
    rgrid->SetYCoordinates(da_x);

    // set point data of structured grid as the vtkDataArray
    vtkPointData * pointData = rgrid->GetPointData();
    pointData->SetScalars(scalars);


    vtkMapper * mapper;

    if (opts.do_contour) {
        if (contour == NULL) 
            contour = vtkContourFilter::New();
        contour->SetInput( rgrid );
        contour->GenerateValues( opts.numContours, min_v, max_v);
        if (polyMapper == NULL)
            polyMapper = vtkPolyDataMapper::New();
        polyMapper->SetInput( contour->GetOutput() );
        mapper = polyMapper;
    } else {
        // Data mapper
        if (dataMapper == NULL) {
            dataMapper = vtkDataSetMapper::New();
        } 
        dataMapper->SetInput(rgrid);
        //dataMapper->SetInput(filter->GetOutput());
        mapper = dataMapper;
    }


    // debug
    if (verbose>1) { 
    printf("    vtk_plotRectlinearGrid:\n");
    printf("      v dims:           %dx%d = %d points\n", dims[0], dims[1], dims[0]*dims[1]);
    printf("      scalars:\n");
    printf("      pointdata:        tuples = %d\n", pointData->GetAbstractArray(0)->GetNumberOfTuples());
    printf("      rgrid:   points = %d\n", rgrid->GetNumberOfPoints());
    }

    scalars->Delete();
    rgrid->Delete();

    return vtk_plot2D( mapper, NULL, opts, dims[0], dims[1], min_v, max_v, min_x, max_x, min_y, max_y);
}


/**
 *   Plot a 2D array over a STRUCTURED GRID
 *   2D x and y arrays should be supplied for a 2D array (with the same dimensions as v).
 */
int vtk_plotStructuredGrid( Array v, double min_v, double max_v, 
                            Array x, double min_x, double max_x,
                            Array y, double min_y, double max_y,
                            vtk_options opts)
{
    int dims[2];
    for (int i=0; i<2; i++) dims[i] = Array_getDimensionSize(v, i);

    // make a vtkDataArray from the Array 
    vtkDataArray * scalars = arrayToVTKArray(v);
    //vtkDataArray * da_x = arrayToVTKArray(x);
    //vtkDataArray * da_y = arrayToVTKArray(y);

    // build the 3D points from x, y and 0
    vtkPoints * points = vtkPoints::New();
    points->SetNumberOfPoints(dims[1]*dims[0]);
    int offset; 
    double dx, dy;
    // v,x,y values are stored row-major order
    // but Points should be stored column-major order
    for (int j=0; j<dims[0]; j++) {
        offset = j*dims[1];
        for (int i=0; i<dims[1]; i++) {
            //dx = da_x->GetTuple1(i+offset);
            //dy = da_y->GetTuple1(i+offset);
            dx = Array_doublevalue(x,i+offset);
            dy = Array_doublevalue(y,i+offset);
            points->SetPoint(i+offset, dx, dy, 0.0);
            //printf(" SetPoint( %d, %f, %f, 0.0)\n", i+offset, dx, dy);
        }
    }
    printf(".... set %dx%d = %d points\n", dims[1], dims[0], dims[1]*dims[0]);
    
    // Structured Grid 
    vtkStructuredGrid * sgrid 
            = vtkStructuredGrid::New();
    sgrid->SetDimensions(dims[1],dims[0],1);
    sgrid->SetPoints(points);
    //sgrid->SetScalarType( scalars->GetDataType() );

    // set point data of structured grid as the vtkDataArray
    vtkPointData * pointData = sgrid->GetPointData();
    pointData->SetScalars(scalars);

    /*
    vtkStructuredGridGeometryFilter * filter =
        vtkStructuredGridGeometryFilter::New();
    filter->SetInput(sgrid);
    filter->SetExtent(0, dims[1], 0, dims[0], 0, 0);
    */

    vtkMapper * mapper;

    if (opts.do_contour) {
        if (contour == NULL) 
            contour = vtkContourFilter::New();
        contour->SetInput( sgrid );
        contour->GenerateValues( opts.numContours, min_v, max_v);
        if (polyMapper == NULL)
            polyMapper = vtkPolyDataMapper::New();
        polyMapper->SetInput( contour->GetOutput() );
        mapper = polyMapper;
    } else {
        // Data mapper
        if (dataMapper == NULL) {
            dataMapper = vtkDataSetMapper::New();
        } 
        dataMapper->SetInput(sgrid);
        //dataMapper->SetInput(filter->GetOutput());
        mapper = dataMapper;
    }


    // debug
    if (verbose>1) { 
    printf("    vtk_plotStructuredGrid:\n");
    printf("      v dims:           %dx%d = %d points\n", dims[0], dims[1], dims[0]*dims[1]);
    printf("      scalars:\n");
    printf("      pointdata:        tuples = %d\n", pointData->GetAbstractArray(0)->GetNumberOfTuples());
    printf("      sgrid:   points = %d\n", sgrid->GetNumberOfPoints());
    }

    points->Delete();
    scalars->Delete();
    sgrid->Delete();

    return vtk_plot2D( mapper, NULL, opts, dims[0], dims[1], min_v, max_v, min_x, max_x, min_y, max_y);
}

/**
 *   Plot a 1D array over a STRUCTURED GRID
 *   1D x and y arrays should be supplied for a 1D array (with the same dimensions as v).
 */
int vtk_plotStructuredGrid1D( Array v, double min_v, double max_v, 
                            Array x, double min_x, double max_x,
                            Array y, double min_y, double max_y,
                            vtk_options opts)
{
    int size;
    size = Array_getDimensionSize(v, 0);

    // make a vtkDataArray from the Array 
    vtkDataArray * scalars = arrayToVTKArray(v);
    //vtkDataArray * da_x = arrayToVTKArray(x);
    //vtkDataArray * da_y = arrayToVTKArray(y);

    // build the 3D points from x, y and 0
    vtkPoints * points = vtkPoints::New();
    points->SetNumberOfPoints(size);
    double dx, dy;
    // v,x,y values are stored row-major order
    // but Points should be stored column-major order
    for (int j=0; j<size; j++) {
            //dx = da_x->GetTuple1(j);
            //dy = da_y->GetTuple1(j);
            dx = Array_doublevalue(x,j);
            dy = Array_doublevalue(y,j);
            points->SetPoint(j, dx, dy, 0.0);
            //printf(" SetPoint( %d, %f, %f, 0.0)\n", j, dx, dy);
    }
    printf(".... set %d points\n", size);
    
    // Structured Grid 
    vtkStructuredGrid * sgrid 
            = vtkStructuredGrid::New();
    sgrid->SetDimensions(size,1,1);
    sgrid->SetPoints(points);
    //sgrid->SetScalarType( scalars->GetDataType() );

    // set point data of structured grid as the vtkDataArray
    vtkPointData * pointData = sgrid->GetPointData();
    pointData->SetScalars(scalars);

    /*
    vtkStructuredGridGeometryFilter * filter =
        vtkStructuredGridGeometryFilter::New();
    filter->SetInput(sgrid);
    filter->SetExtent(0, dims[1], 0, dims[0], 0, 0);
    */

    vtkMapper * mapper;

    if (opts.do_contour) {
        if (contour == NULL) 
            contour = vtkContourFilter::New();
        contour->SetInput( sgrid );
        contour->GenerateValues( opts.numContours, min_v, max_v);
        if (polyMapper == NULL)
            polyMapper = vtkPolyDataMapper::New();
        polyMapper->SetInput( contour->GetOutput() );
        mapper = polyMapper;
    } else {
        // Data mapper
        if (dataMapper == NULL) {
            dataMapper = vtkDataSetMapper::New();
        } 
        dataMapper->SetInput(sgrid);
        //dataMapper->SetInput(filter->GetOutput());
        mapper = dataMapper;
    }


    // debug
    if (verbose>1) { 
    printf("    vtk_plotStructuredGrid1D:\n");
    printf("      v dims:           %d points\n", size);
    printf("      scalars:\n");
    printf("      pointdata:        tuples = %d\n", pointData->GetAbstractArray(0)->GetNumberOfTuples());
    printf("      sgrid:   points = %d\n", sgrid->GetNumberOfPoints());
    }

    points->Delete();
    scalars->Delete();
    sgrid->Delete();

    return vtk_plot2D( mapper, NULL, opts, size, 1, min_v, max_v, min_x, max_x, min_y, max_y);
}



/**
 *   Plot a 2D array with POLAR COORDINATES (uses Structured Grid).
 *   1D r and phi arrays should be supplied.
 */

int vtk_plotPolarGrid( Array v,   double min_v,   double max_v, 
                       Array r,   double min_r,   double max_r,
                       Array phi, double min_phi, double max_phi,
                       vtk_options opts)
{
    int dims[2];
    /*
    Array r, phi;
    if (r_phi_pos == 1) {
        r   = x;
        phi = y;
    } else {
        r   = y;
        phi = x;
    }
    */
    int rdim = Array_getDimensionSize(r,0);
    int phidim = Array_getDimensionSize(phi,0);
    for (int i=0; i<2; i++) dims[i] = Array_getDimensionSize(v, i);

    // make a vtkDataArray from the Array 
    vtkDataArray * scalars = arrayToVTKArray(v);

    // build the 3D points from r, phi and 0
    vtkPoints * points = vtkPoints::New();
    points->SetNumberOfPoints(rdim*phidim);
    int offset; 
    double rv, phiv;
    if (opts.radius_is_first_dim == 1) {
        for (int j=0; j<rdim; j++) { // 
            rv = Array_doublevalue(r,j);
            offset = j*phidim;
            for (int i=0; i<phidim; i++) {
                phiv = Array_doublevalue(phi,i);
                points->SetPoint(i+offset, rv*cos(phiv), rv*sin(phiv), 0.0);
                //printf("%5d: r=%f  phi=%f  x=%f  y=%f\n",i+offset, rv, phiv, rv*cos(phiv), rv*sin(phiv), 0.0);
            }
        }
    } else {
        for (int j=0; j<phidim; j++) { // 
            phiv = Array_doublevalue(phi,j);
            offset = j*rdim;
            for (int i=0; i<rdim; i++) {
                rv = Array_doublevalue(r,i);
                points->SetPoint(i+offset, rv*cos(phiv), rv*sin(phiv), 0.0);
                //printf("%5d: r=%.15le  phi=%.15le  x=%f  y=%f\n",i+offset, rv, phiv, rv*cos(phiv), rv*sin(phiv), 0.0);
            }
        }
    }
    if (verbose>2) printf(".... set %dx%d = %d points\n", dims[1], dims[0], dims[1]*dims[0]);
    
    // Structured Grid 
    vtkStructuredGrid * sgrid 
            = vtkStructuredGrid::New();
    sgrid->SetDimensions(dims[1],dims[0],1);
    sgrid->SetPoints(points);

    // set point data of structured grid as the vtkDataArray
    vtkPointData * pointData = sgrid->GetPointData();
    pointData->SetScalars(scalars);

    vtkMapper * mapper;

    if (opts.do_contour) {
        if (contour == NULL) 
            contour = vtkContourFilter::New();
        contour->SetInput( sgrid );
        contour->GenerateValues( opts.numContours, min_v, max_v);
        if (polyMapper == NULL)
            polyMapper = vtkPolyDataMapper::New();
        polyMapper->SetInput( contour->GetOutput() );
        mapper = polyMapper;
    } else {
        // Data mapper
        if (dataMapper == NULL) {
            dataMapper = vtkDataSetMapper::New();
        } 
        dataMapper->SetInput(sgrid);
        mapper = dataMapper;
    }


    // debug
    if (verbose>1) { 
    printf("    vtk_plotPolarGrid\n");
    printf("      v   dims:     %dx%d = %d points\n", dims[0], dims[1], dims[0]*dims[1]);
    printf("      r   dims:     %d\n", Array_getDimensionSize(r,0));
    printf("      phi dims:     %d\n", Array_getDimensionSize(phi,0));
    printf("      scalars:\n");
    printf("      pointdata: tuples = %d\n", pointData->GetAbstractArray(0)->GetNumberOfTuples());
    printf("      sgrid:     points = %d\n", sgrid->GetNumberOfPoints());
    }

    points->Delete();
    scalars->Delete();
    sgrid->Delete();

    return vtk_plot2D( mapper, NULL, opts, dims[0], dims[1], min_v, max_v, -1.0, -1.0, -1.0, -1.0);
}

/**
 *   Plot a 1D array on a TRIANGLE MESH
 *   if length(v) is N, 
 *      coords is a Nx2 or Nx3 array of point x-y[-z] co-ordinates
 *      connections is an int32 or int16 Mx3 array of triangles (with point indexes (to coords array)) 
 */
int vtk_plotTriangleMesh( Array v,           double min_v,           double max_v, 
                          Array coords,      double min_coords,      double max_coords,
                          Array connections, double min_connections, double max_connections,
                          vtk_options opts)
{
    struct timeval tp;
    _tstart();
    int dims[2];
    for (int i=0; i<2; i++) dims[i] = Array_getDimensionSize(v, i);

    // make a vtkDataArray from the value Array v
    vtkDataArray * scalars = arrayToVTKArray(v);

    // build the 3D points of vertices from coords 
    int numpoints = Array_getDimensionSize(coords,0);  // number of points
    int pointdim = Array_getDimensionSize(coords,1);   // 2 or 3 coords of a point is given?
    vtkDoubleArray * pcoords = vtkDoubleArray::New();
    pcoords->SetNumberOfComponents(3);           // 3 components represent a point
    pcoords->SetNumberOfTuples(numpoints);       // allocate space for the points
    double c[3] = { 0.0, 0.0, 0.0};              // coords, c[2] may not be changed at all
    int idx=0;
    for (int i=0; i<numpoints; i++) {
        c[0] = Array_doublevalue(coords, idx++);
        c[1] = Array_doublevalue(coords, idx++);
        if (pointdim == 3) 
            c[2] = Array_doublevalue(coords, idx++);
        pcoords->SetTuple(i, c);
    }
    vtkPoints * points = vtkPoints::New();
    points->SetData(pcoords);

    // build the triangles from connections 
    vtkCellArray* triangles = vtkCellArray::New();
    
    int nTriangles = Array_getDimensionSize(connections, 0);
    vtkIdType p[3];
    int32_t *i32a;
    int16_t *i16a;
    idx = 0;
    if (Array_getType(connections) == int32Array) {
        i32a = (int32_t *) Array_getDataPointer(connections);
        for (int i = 0; i < nTriangles; i++) {
            //idx=3*i;
            p[0] = (vtkIdType) i32a[idx++];
            p[1] = (vtkIdType) i32a[idx++];
            p[2] = (vtkIdType) i32a[idx++];
            triangles->InsertNextCell( 3, p);
        }
    } else if (Array_getType(connections) == int16Array) {
        i16a = (int16_t *) Array_getDataPointer(connections);
        for (int i = 0; i < nTriangles; i++) {
            p[0] = (vtkIdType) i16a[idx++];
            p[1] = (vtkIdType) i16a[idx++];
            p[2] = (vtkIdType) i16a[idx++];
            triangles->InsertNextCell( VTK_TRIANGLE, p);
        }
    } else {
        fprintf(stderr, "Error for triangle mesh: Array y must be a 16 or 32 bit integer array (indexes of x's values)\n");
        return 1;

    }
    

    // Create the triangle mesh dataset. 
    vtkPolyData* polydata = vtkPolyData::New();
    // Assign points and cells
    polydata->SetPoints(points);
    polydata->SetPolys(triangles);
    // Assign scalars
    polydata->GetPointData()->SetScalars(scalars);
        
    // Create the polymapper
    if (opts.do_contour) {
        if (contour == NULL) 
            contour = vtkContourFilter::New();
        contour->SetInput( polydata );
        contour->GenerateValues( opts.numContours, min_v, max_v);
        if (polyMapper == NULL)
            polyMapper = vtkPolyDataMapper::New();
        polyMapper->SetInput( contour->GetOutput() );
    } else {
        if (polyMapper == NULL)
            polyMapper = vtkPolyDataMapper::New();
        polyMapper->SetInput( polydata );
    }


    // debug
    if (verbose>1) { 
    printf("    vtk_plotTriangleMesh\n");
    printf("      value     dims:     %d points\n", Array_getDimensionSize(v, 0));
    printf("      vertices  num:      %d\n", Array_getDimensionSize(coords,0));
    printf("      triangles num:      %d\n", Array_getDimensionSize(connections,0));
    printf("      scalars:\n");
    printf("      pointdata: tuples    = %d\n", polydata->GetPointData()->GetAbstractArray(0)->GetNumberOfTuples());
    printf("      polydata:  points    = %d\n", polydata->GetNumberOfPoints());
    printf("      polydata:  triangles = %d\n", polydata->GetNumberOfPolys());
    }

    points->Delete();
    pcoords->Delete();
    scalars->Delete();
    triangles->Delete();
    polydata->Delete();
    
    if (verbose>1) {
        _tend();
        _tdiff(&tp);
        printf("     time in vtk_plotTriangleMesh = %d.%6.6d sec\n", tp.tv_sec, tp.tv_usec);
    }
    return vtk_plot2D(polyMapper, NULL, opts, -1, -1, min_v, max_v, -1.0, -1.0, -1.0, -1.0);
}

/** Make a 2D plot offscreen and write into a file. 
 *  Mapper should be the input to the pipeline. 
 *  dim_x/dim_y can be -1,-1 if it is meaningless for a plot (e.g. triangle mesh)
 *  min_x/max_x can be -1,-1 if it is meaningless for a plot (e.g. triangle mesh)
 *  min_y/max_y can be -1,-1 if it is meaningless for a plot (e.g. triangle mesh)
 *
 *  On the image, X is the vertical scale and Y is the horizontal, because C arrays
 *  are stored row-major order while VTK represent points in column-major order.
 */
static int vtk_plot2D( vtkMapper * mapper, vtkCubeAxesActor2D * axesActor, vtk_options opts, 
                       int dim_x,    int dim_y, 
                       double min_v, double max_v,
                       double min_x, double max_x,
                       double min_y, double max_y)
{
    struct timeval tp;
    _tstart();
    // Get selected colormap
    vtkLookupTable * lut = getSelectedLookupTable(min_v, max_v);

#define ABS(x) (x >= 0.0 ? x : -1.0*x)
    if (opts.balance_colormap && min_v < 0.0 && max_v > 0.0) {
        double dist = max_v - min_v;
        double diff = max_v + min_v;
        diff = ABS(diff);
        if ( diff/dist < 0.5 ) {
            if      (min_v > -1.0*max_v) min_v = -1.0*max_v;
            else if (max_v < -1.0*min_v) max_v = -1.0*min_v;
            if (verbose>1) printf("      balanced colormap to %lf..%lf\n", min_v, max_v);
        } else {
            printf("Warning: Colormap will not be balanced for %s because difference of"
                    "min value =%lf and max value = %lf is more then 50%% of their distance = %lf\n",
                    opts.title, min_v, max_v, dist);
        }
    }
    mapper->SetScalarRange(min_v, max_v);
    mapper->SetLookupTable(lut);


    // Actor
    if (verbose>2) printf("      make actor...\n");
    if (dataActor == NULL) {
        dataActor = vtkActor::New();
    }
    dataActor->SetMapper(mapper);

    // scale both sides of the array view to the same size (make square image)
    if (opts.do_square) {
        double hscale = 1.0, vscale = 1.0;
        if (dim_x < dim_y) 
            vscale = double(dim_y) / double(dim_x);
        else
            hscale = double(dim_x) / double(dim_y);
        dataActor->SetScale( hscale, vscale, 1);
    }
    
    // create color legend
    if (scalarBar == NULL)  {
        scalarBar = vtkScalarBarActor::New();
        // configure the legend font
        scalarBar->GetLabelTextProperty()->SetBold(0);
        scalarBar->GetLabelTextProperty()->SetShadow(0);
        scalarBar->GetLabelTextProperty()->SetColor(foregroundRGB);
        scalarBar->GetTitleTextProperty()->SetBold(0);
        scalarBar->GetTitleTextProperty()->SetItalic(0);
        scalarBar->GetTitleTextProperty()->SetShadow(0);
        scalarBar->GetTitleTextProperty()->SetColor(foregroundRGB);
        //scalarBar->GetProperty()->SetColor(0,0,0);
        //scalarBar->GetProperty()->SetOpacity(0.0);

        // configure position, orientation and size
        scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        scalarBar->GetPositionCoordinate()->SetValue( 0.1, 0.02 );
        scalarBar->SetOrientationToHorizontal();
        scalarBar->SetTextPositionToPrecedeScalarBar();
        scalarBar->SetWidth(0.8);
        scalarBar->SetHeight(0.98);
    }
    scalarBar->SetLookupTable( mapper->GetLookupTable() );
    scalarBar->SetTitle(opts.title);

    // create axes
    /*
    if (opts.show_axes) {
        double b[6];
        dataActor->GetBounds(b);
        //double bounds[6]={min_x, max_x, min_y, max_y, 0.0, 0.0};
        //double range[6]={min_v, max_v, min_v, max_v, 0.0, 0.0};
        // Note that axes are flipped for the arrays (2D C array -> VTK 2D array)
        double range[6]={min_y, max_y, min_x, max_x, 0.0, 0.0};
        printf("### Axis bounds: x=%f..%f  y=%f..%f\n", b[0], b[1], b[2], b[3]);
        printf("### Axis ranges: x=%f..%f  y=%f..%f\n", min_y, max_y, min_x, max_x);
        if (axesActor == NULL) {
            axesActor = vtkCubeAxesActor2D::New();
            axesActor->SetBounds(b);
            axesActor->SetRanges(range);
            axesActor->SetUseRanges(1);
            axesActor->SetZAxisVisibility(0);
            vtkTextProperty * txtp = vtkTextProperty::New();
            txtp->SetColor(0,0,0);
            axesActor->SetAxisLabelTextProperty(txtp);
            axesActor->SetAxisTitleTextProperty(txtp);
            axesActor->SetFontFactor(1.4); // txtp->SetFontSize has no effect 
            axesActor->SetNumberOfLabels(4);
            //printf("label format = %s\n", axesActor->GetLabelFormat());
            axesActor->SetLabelFormat("%.0f");
        }
        // Note that axis names are flipped too
        if (opts.xname != NULL) 
            axesActor->SetYLabel(opts.xname);
        else
            axesActor->SetYLabel("x");
        if (opts.yname != NULL) 
            axesActor->SetXLabel(opts.yname);
        else
            axesActor->SetXLabel("y");

        
    }
    */

#ifndef POPUP_RENDERWINDOW
    // ensure drawing offscreen
    vtk_set_offscreen();
#endif
    
    if (verbose>2) printf("      create renderers...\n");
    // Renderers
    if (renderer == NULL) {
        renderer = vtkRenderer::New();
        renderer->SetBackground( backgroundRGB );
        renderer->SetViewport( 0.0, 0.1, 1.0, 1.0);
        renderer->AddActor( dataActor );
    }
    
    renderer->ResetCamera();
    vtkCamera * camera = renderer->GetActiveCamera();
    if (axesActor != NULL) {
        double b[6];
        dataActor->GetBounds(b);
        axesActor->SetBounds(b);
        //printf("### Axis bounds: x=%f..%f  y=%f..%f\n", b[0], b[1], b[2], b[3]);
        axesActor->SetCamera(camera);
        renderer->AddActor(axesActor);

    }
    
    if (!opts.show_axes && opts.auto_zoom) {
        double scale = 1.0;
        double hlen = max_y - min_y;  // horizontal (width) is given in y
        double vlen = max_x - min_x;  // vertical (height) is given in x
        if (opts.do_square || hlen == vlen) 
            scale = 1.0; // if we scaled the actor to square plot, zoom as for x=y
        else if (hlen < vlen) 
            scale = double(vlen) / double(hlen);
        else
            scale = double(hlen) / double(vlen);
        // zoom: 1.3 at x=y (square, scale=1.0), 1.0 at scale=2.0, and do not go below 1.0
        scale = (scale < 2.0 ? scale : 2.0);
        camera->Zoom(1.3 - 0.3*(scale-1.0));
        if (verbose>2) printf("      zoom: scale = %lf, zoom = %lf\n", scale, 1.3 - 0.3*(scale-1.0));
    }
    
    
    if (barrenderer == NULL) {
        barrenderer = vtkRenderer::New();
        barrenderer->SetBackground( backgroundRGB );
        barrenderer->SetViewport( 0.0, 0.0, 1.0, 0.1);
        barrenderer->AddActor2D( scalarBar );
    }
    
    if (verbose>2) printf("      create render window...\n");
    // Render Window
    if (render_window == NULL) {
        render_window = vtkRenderWindow::New();
#ifndef POPUP_RENDERWINDOW
        render_window->SetOffScreenRendering( 1 );
#endif
        render_window->AddRenderer( renderer );
        render_window->AddRenderer( barrenderer );
    }
    render_window->SetSize( imgSizeX, imgSizeY );

    // Window to Image        
    if (win_2_image == NULL) { 
        win_2_image = vtkWindowToImageFilter::New();
        win_2_image->SetInput( render_window );
    }

    if (verbose>2) printf("      create writer...\n");
    // PNG Writer
    char outfile[256];
    snprintf(outfile, 256, "%s.png", opts.imagefilename);
    if (pngwriter == NULL) {
        pngwriter = vtkPNGWriter::New();
        pngwriter->SetInputConnection( win_2_image->GetOutputPort() );
    }
    pngwriter->SetFileName( outfile );

    // Do rendering and write to file
    if (verbose>2) printf("      do rendering...\n");
    if (verbose>1) {
        _tend();
        _tdiff(&tp);
        printf("     time in vtk_plot2D section 2 = %d.%6.6d sec\n", tp.tv_sec, tp.tv_usec);

        _tstart();
    }
    render_window->Render();

    if (verbose>2) printf("      do writing...\n");
    win_2_image->Modified();
    if (verbose>1) {
        _tend();
        _tdiff(&tp);
        printf("     time in vtk_plot2D section 3 = %d.%6.6d sec\n", tp.tv_sec, tp.tv_usec);

        _tstart();
    }
    pngwriter->Write();
    if (verbose>1) {
        _tend();
        _tdiff(&tp);
        printf("     time in vtk_plot2D/pngwriter->Write() = %d.%6.6d sec\n", tp.tv_sec, tp.tv_usec);
    }

#ifdef POPUP_RENDERWINDOW
    sleep(1.0);
#endif

    if (print_provenance)
        printf("--ProvenanceInfo name=\"%s\" size=0 date=0\n", outfile);

    return 0;
}



