#ARCH=macosx
#ARCH=ewok
#ARCH=xt5
ARCH=sith
#ARCH=titan

# OPTIONS TO BE SET SHOULD BE ABOVE THIS LINE

ifdef RECURSIVE_CALL
include ../Makefile.${ARCH}
else
include Makefile.${ARCH}
endif

ifeq ("${USE_VTK}","yes")
PROG=plotter plotter2d
else
PROG=plotter
endif

default: ${PROG}
all:     ${PROG}

plotter: build.plotter array.c dirutil.c minmax.c reader.c main.c ${SRC_NETCDF} ${SRC_HDF5} ${SRC_ADIOS} ${SRC_XMGRACE}
	@echo " "
	@echo "----- PLOTTER ------"
	(cd build.plotter; RECURSIVE_CALL=yes make -f ../Makefile plotter_exe)

plotter2d: build.plotter2d array.c dirutil.c minmax.c reader.c main.c ${SRC_NETCDF} ${SRC_HDF5} ${SRC_ADIOS} ${SRC_XMGRACE} ${SRC_VTK}
	@echo " "
	@echo "----- PLOTTER 2D ------"
	(cd build.plotter2d; RECURSIVE_CALL=yes \
	ADD_CFLAGS="$(VTK_CFLAGS)" \
	ADD_CXXFLAGS="$(VTK_CXXFLAGS)" \
	ADD_INC="$(VTK_INC)" \
	ADD_LIBS="$(VTK_LIBS)" \
	ADD_LDFLAGS="$(VTK_LDFLAGS)" \
	make -f ../Makefile plotter2d_exe)

build.plotter:
	mkdir -p build.plotter

build.plotter2d:
	mkdir -p build.plotter2d

# THE RULES BELOW WILL BE USED IN SUBDIRS build.plotter*

# $(<D) = The directory part of the first prerequisite.
# $(<F) = The file-within-directory part of the first prerequisite

%.o: ../%.c
	$(CC) $(CFLAGS) $(ADD_CFLAGS) -c  ${INC} ${ADD_INC} $<

%.o: ../%.cxx
	$(CXX) $(CXXFLAGS) $(ADD_CXXFLAGS) -c  ${INC} ${ADD_INC} $< 

plotter_exe: array.o dirutil.o minmax.o reader.o ${OBJ_NETCDF} ${OBJ_HDF5} ${OBJ_ADIOS} ${OBJ_XMGRACE} timer.o main.o
	$(CC) $(LDFLAGS) -o ../plotter array.o dirutil.o minmax.o reader.o ${OBJ_NETCDF} ${OBJ_HDF5} ${OBJ_ADIOS} ${OBJ_XMGRACE} timer.o main.o ${LIBS}

plotter2d_exe: array.o dirutil.o minmax.o reader.o ${OBJ_NETCDF} ${OBJ_HDF5} ${OBJ_ADIOS} ${OBJ_XMGRACE} ${OBJ_VTK} timer.o main.o
	$(CXX) $(LDFLAGS) $(ADD_LDFLAGS) -o ../plotter2d array.o dirutil.o minmax.o reader.o ${OBJ_NETCDF} ${OBJ_HDF5} ${OBJ_ADIOS} ${OBJ_XMGRACE} ${OBJ_VTK} timer.o main.o ${LIBS} ${ADD_LIBS}

clean:
	rm -f ${PROG} *.o *~
	rm -f testArray 
	rm -rf build.plotter build.plotter2d

h5_read: h5_read.o
	$(CC) $(CFLAGS) -o h5_read h5_read.o ${LIBS}

h5trav: h5trav.o
	$(CC) $(CFLAGS) -o h5trav h5trav.o ${LIBS}

