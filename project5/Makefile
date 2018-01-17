################################################################################
#                                                                              #
#                                   MAKEFILE                                   #
#                                                                              #
#                                                                              # 
#     COMPILE ALL TARGETS: ...... "make" (without the "")                      #
#     COMPILE <program.cc>: ..... "make program"                               #
#     CLEAN ALL EXEC/OBJ-FILES:.. "clean"                                      #
#                                                                              #
################################################################################

# CCFLAGS: -W-Flags simply control the warning output
#          -O3 is rougly speaking for speeding up the code - of course, this makes the compilation process a teeny tiny tad bit slower
#          -std=c++11 defines the compiler standard to adhere to c++11, which contains many additional features as compared to c++
CCFLAGS = -Wall -Wextra -Wno-deprecated -Wno-reorder -Wno-unused-local-typedefs -O3 -std=c++11
CC = /usr/lib64/mpich-3.2/bin/mpicc
# INCFLAGS: Many libraries include files in the header without specifying the exact location, despite the .h files not being in the same directory or one of the standard system directories. In this case, one has to tell the compiler, which directories to search through, which is accomplished using -I
INCFLAGS = -I$(GEN_INCLUDE) -I$(BOOST_INCLUDE) -I$(EIGEN_DIR) -I$(MPI_INCLUDE) -I$(PETSC_INCLUDE) -I$(SRC_DIR) -I$(STLIB_OLD_INCLUDE) -I$(SUPERLU_INCLUDE) -I$(VTK_INCLUDE)
# -I$(STLIB_INCLUDE) 

# LDFLAGS: Oftentimes, libraries are not included as .h 'human-readable-files', but as .o or .so's binary files (e.g. to keep the code secret). In this case, we ask to bind in a binary library via -lsomename, which would look for a file libsomename.so
LDFLAGS = -lblas -lpetsc -lelemental -lpmrrr -lsuperlu_dist_3.3 -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lX11 -lparmetis -lmetis -lpthread -lmpichf90 -lgfortran -lm -lquadmath -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lmpl -lrt -lgcc_s -lXt -lGLU -lGL -lSM -lICE -lXext

# LIBFLAGS Oftentimes, the libraries called in LDFLAGS are not present in the current working directory, nor in one of the standard system directories, in which case, we need to tell the compiler where to look. This is accomplished using the -L flag.
LIBFLAGS = -L$(GEN_LIB) -L$(GEN_LIB64) -L$(MPI_LIB) -L$(PETSC_LIB) -L$(VTK_LIB)

# RENDERINGFLAGS Just some rendering options, which I'm not even sure if they're of any use, really...
RENDERINGFLAGS = -DvtkRenderingCore_AUTOINIT="4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)" -DvtkRenderingVolume_AUTOINIT="1(vtkRenderingVolumeOpenGL)"

# RPATH In some cases, libraries are not only required at compilation time (i.e. the stuff that's happening when you write 'make' or 'make Main') but at runtime, in which we need to specify the runtime search path (hence the designation RPATH)
RPATH = -Wl,-rpath,$(MPI_LIB) -Wl,-rpath,$(PETSC_LIB) -Wl,-rpath,$(VTK_LIB)  

# HEADER FILES Simply defines which '.h''s to "possibly" include
HEADER_FILES = $(wildcard *.h)

# TARGETS File(s) to be compiled
TARGETS = Main_Task3 Main_Task4

# NOTE The following lines are pretty standard and the meaning of shortcuts such as @ (don't echo - check it yourself...if you "make", you won't see any message apart from the one defined in the line above), $@, $<, ${...}, all, clean, distclean, again, and so on can be found here: http://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html
  
all: $(TARGETS)

%: %.cc $(wildcard *.h)
	@echo "Building $@"
	@${CC} ${CCFLAGS} -o $@ $< ${INCFLAGS} ${LDFLAGS} ${LIBFLAGS} ${RENDERINGFLAGS} ${RPATH}  ${VTK_LIB}/*.so.*

clean:
	$(RM) $(TARGETS) *~ core* *.stackdump *.o

distclean:
	$(MAKE) clean
	$(RM) $(TARGETS)

again:
	$(MAKE) distclean
	$(MAKE)
