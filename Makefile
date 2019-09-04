# -----------------------------------------------------------------
# Version: 4.0
# Date: July 2018 
# Makefile for PIHM v4.0
# -----------------------------------------------------------------
# Programmer: Lele Shu (lele.shu@gmail.com)
# -----------------------------------------------------------------
#  Prerequisite:
#  1 install sundials 3.0+ via https://computation.llnl.gov/projects/sundials/sundials-software.
#  2 If parallel-computing is prefered, please install OpenMP.
#	 For mac: 
#	  		brew install llvm clang
#			brew install libomp
#			compile flags for OpenMP: 
#				-Xpreprocessor -fopenmp -lomp
#			Library/Include paths:
#				-L/usr/local/opt/libomp/lib 
#				-I/usr/local/opt/libomp/include
#			
# -----------------------------------------------------------------
# Configure this File:
# 1 Path of SUNDIALS_DIR. [CRITICAL]
# 2 Path of OpenMP if parallel is preffered.
# 3 Path of SRC_DIR, default is "SRC_DIR = ."
# 4 Path of BUILT_DIR, default is "BUILT_DIR = ."
# -----------------------------------------------------------------
SUNDIALS_DIR = /usr/local/sundials
# SUNDIALS_DIR = ~/sundials


SHELL = /bin/sh
BUILDDIR = .
SRC_DIR = src

LIB_SYS = /usr/local/lib/
INC_OMP = /usr/local/opt/libomp/include
LIB_OMP = /usr/local/opt/libomp/lib

INC_MPI = /usr/local/opt/open-mpi

TARGET_PIHM     = ${BUILDDIR}/pihm++
TARGET_OMP      = ${BUILDDIR}/pihm_omp
TARGET_CB_MPI   = ${BUILDDIR}/calib_mpi
TARGET_CB_OMP   = ${BUILDDIR}/calib_omp
TARGET_DEBUG    = ${BUILDDIR}/pihm_debug

MAIN_PIHM 		= ${SRC_DIR}/PIHMmain.cpp
MAIN_OMP 		= ${SRC_DIR}/PIHMmain.cpp
MAIN_CB_MPI 	= ${SRC_DIR}/PIHMcalib_MPI.cpp
MAIN_CB_OMP 	= ${SRC_DIR}/PIHMcalib_OpenMP.cpp
MAIN_DEBUG 		= ${SRC_DIR}/PIHMmain.cpp

# If compile on Cluster
# CC       = g++
# MPICC    = mpic++
# LK_OMP   = -fopenmp -lsundials_nvecopenmp
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SUNDIALS_DIR}/lib

CC       = /usr/bin/g++
MPICC    = /usr/local/bin/mpic++
CFLAGS   = -O3 -g  -std=c++11
#STCFLAG     = -static
LDFLAGS  = 
LIBS     = -lm
SRC    	= ${SRC_DIR}/classes/*.cpp \
		  ${SRC_DIR}/ModelData/*.cpp \
		  ${SRC_DIR}/Model/*.cpp \
		  ${SRC_DIR}/Equations/*.cpp

SRC_H	= ${SRC_DIR}/classes/*.hpp \
		  ${SRC_DIR}/ModelData/*.hpp \
		  ${SRC_DIR}/Model/*.hpp \
		  ${SRC_DIR}/Equations/*.hpp


INCLUDES = -I ${SUNDIALS_DIR}/include \
		   -I ${SUNDIALS_DIR}/include/cvode \
		   -I ${SUNDIALS_DIR}/include/nvector \
		   -I ${SUNDIALS_DIR}/include/sundials \
		   -I ${INC_OMP} \
		   -I ${SRC_DIR}/Model \
		   -I ${SRC_DIR}/ModelData \
		   -I ${SRC_DIR}/classes \
		   -I ${SRC_DIR}/Equations 

SRC_CB = ${SRC} \
			${SRC_DIR}/PIHMcalib/*.cpp
SRC_CB_H = ${SRC_H} \
			  ${SRC_DIR}/PIHMcalib/*.hpp
INC_CB = ${INCLUDES} \
		   -I ${SRC_DIR}/PIHMcalib
		  
LIBRARIES = -L ${LIB_OMP} \
			-L ${SUNDIALS_DIR}/lib \
			-L ${LIB_SYS}

LK_FLAGS = -lm -lsundials_cvode -lsundials_nvecserial
LK_OMP	= -Xpreprocessor -fopenmp -lomp -lsundials_nvecopenmp
LK_CALIB = _CALIBMODE



all:
	@echo 
	@echo 'make pihm'
	make pihm
	@echo
	
	@echo 'make pihm_omp'
	make pihm_omp
	@echo
	@echo 'make calib_mpi'
	make calib_mpi
	@echo
	
	@echo 'make calib_omp'
	make calib_omp
	@echo
	
help:
	@(echo)
	@echo "Usage:"
	@(echo '       make all	    	- make both pihm and pihm_omp')
	@(echo '       make pihm     	- make pihm executable')
	@(echo '       make pihm_omp    - make pihm_omp with OpenMP support')
	@(echo '       make calib_mpi   - make calib_mpi with OpenMPI support')
	@(echo '       make calib_omp   - make calib_omp with OpenMP support')
	@(echo)
	@(echo '       make clean    	- remove all executable files')
	@(echo)

pihm: ${MAIN_PIHM} $(SRC) $(SRC_H)
	@echo '...Compiling PIHM ...'
	@echo $(CC) $(CFLAGS) ${STCFLAG} ${INCLUDES} ${LIBRARIES} -o ${TARGET_PIHM} ${MAIN_PIHM} $(SRC)  $(LK_FLAGS)
	@echo
	@echo
	$(CC) $(CFLAGS) ${INCLUDES} ${STCFLAG} ${LIBRARIES} -o ${TARGET_PIHM} ${MAIN_PIHM} $(SRC)  $(LK_FLAGS)
	@echo
	@echo
	@echo " ${TARGET_PIHM} is compiled successfully!"
	@echo

pihm_omp: ${MAIN_OMP}  $(SRC) $(SRC_H)
	@echo '...Compiling PIHM_OpenMP ...'
	@echo $(CC) $(CFLAGS) ${STCFLAG} -D_PIHMOMP ${INCLUDES} ${LIBRARIES} -o ${TARGET_OMP}   ${MAIN_OMP} $(SRC)  $(LK_FLAGS) $(LK_OMP)
	@echo
	@echo
	$(CC) $(CFLAGS)  ${STCFLAG} -D_PIHMOMP ${INCLUDES} ${LIBRARIES} -o ${TARGET_OMP}   ${MAIN_OMP} $(SRC)  $(LK_FLAGS) $(LK_OMP)
	@echo
	@echo " ${TARGET_OMP} is compiled successfully!"
	@echo
	@echo

calib_mpi: ${MAIN_CB_MPI}  $(SRC_CB) $(SRC_CB_H)
	@echo '...Compiling PIHM_Calib OpenMPI...'
	@echo $(MPICC) $(CFLAGS) -D_CALIBMODE ${INC_CB} -I ${INC_MPI} ${LIBRARIES} -o ${TARGET_CB_MPI} ${MAIN_CB_MPI}  $(SRC_CB)  $(LK_FLAGS)
	@echo
	@echo
	 $(MPICC) $(CFLAGS) -D_CALIBMODE ${INC_CB} -I ${INC_MPI} ${LIBRARIES} -o ${TARGET_CB_MPI} ${MAIN_CB_MPI}  $(SRC_CB)  $(LK_FLAGS)
	@echo
	@echo " ${TARGET_CB_MPI} is compiled successfully!"
	@echo
	@echo
calib_omp: ${MAIN_CB_OMP}  $(SRC_CB) $(SRC_CB_H)
	@echo '...Compiling PIHM_Calib OpenMP...'
	@echo $(CC) $(CFLAGS) ${STCFLAG} -D_CALIBMODE ${INC_CB} -I ${INC_MPI} ${LIBRARIES} -o ${TARGET_CB_OMP} ${MAIN_CB_OMP}  $(SRC_CB)  $(LK_FLAGS) $(LK_OMP)
	@echo
	@echo
	$(CC) $(CFLAGS) ${STCFLAG} -D_CALIBMODE ${INC_CB} -I ${INC_MPI} ${LIBRARIES} -o ${TARGET_CB_OMP} ${MAIN_CB_OMP}  $(SRC_CB)  $(LK_FLAGS) $(LK_OMP)
	@echo
	@echo " ${TARGET_CB_OMP} is compiled successfully!"
	@echo
	@echo


clean:
	@echo "Cleaning ... "
	@echo
	@echo "  rm -f *.o"
	@rm -f *.o
	
	@echo "  rm -f ${TARGET_PIHM}"
	@rm -f ${TARGET_PIHM}
	
	@echo "  rm -f ${TARGET_OMP}"
	@rm -f ${TARGET_OMP}
	
	@echo "  rm -f ${TARGET_CB_MPI}"
	@rm -f ${TARGET_CB_MPI}
	
	@echo "  rm -f ${TARGET_CB_OMP}"
	@rm -f ${TARGET_CB_OMP}
	
	@echo
	@echo "Done."
	@echo





