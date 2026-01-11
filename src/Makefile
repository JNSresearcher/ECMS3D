F90=gfortran
CC =gcc

# F90=ifx 
# F90=ifort

SP=./SPARSKIT2
F=ECMS3D 
DLL="C:\Program Files (x86)\Intel\oneAPI\mkl\2023.2.0\redist\intel64"

ifeq ($(F90),gfortran) 

TARGET =$(F)
# -fbounds-check  -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized 
FCOPT1 = -O2  -finit-local-zero -ffree-line-length-none   
FCOPT2 = 
FCOPT3 = -O2 
# -cpp
Lsw = -lz -lgfortran -lstdc++  $(DLL)\mkl_rt.2.dll  $(DLL)\mkl_core.2.dll $(DLL)\mkl_intel_thread.2.dll 
SRC =  uncompress  m_vxc2data m_fparser vxc2data $(SP)\ilut solvers utilites  $(F)
OBJ = $(addsuffix .o, $(SRC) )

.SUFFIXES: .o .mod .f90 .F90 .f .cpp
%.o: %.f90
	$(F90) -o $@ -c $< $(FCOPT1) 
%.o: %.f
	$(F90) -o $@ -c $< $(FCOPT3) 
%.o: %.cpp
	$(CC)  -o $@ -c $< $(FCOPT2) 
$(TARGET): $(OBJ)
	$(F90)  -o $@  $^  $(Lsw)
endif

ifeq ($(F90),ifx) 
# /check:bounds 
TARGET =$(F)
FCOPT1 = /QO2  /Qmkl /heap-arrays:2000 /assume:buffered_io -nologo 
Lsw =  /Qmkl /heap-arrays:2000 
SRC = m_vxc2data m_fparser vxc2data_ifort solvers $(SP)\ilut utilites $(F)
OBJ = $(addsuffix .obj, $(SRC) )
.SUFFIXES: .obj .mod .f90 .F90 .f
%.obj: %.f90
	$(F90) -o $@ -c $< $(FCOPT1)  
%.obj: %.f
	$(F90) -o $@ -c $< $(FCOPT3) 
$(TARGET): $(OBJ)
	$(F90)  -o $@  $^  $(Lsw)
endif

ifeq ($(F90),ifort) 
# /check:bounds 
TARGET =$(F)
FCOPT1 = /QO2  /Qzero /Qmkl /heap-arrays:2000 /assume:buffered_io -nologo 
Lsw =  /Qmkl /heap-arrays:2000 
SRC = m_vxc2data m_fparser vxc2data_ifort solvers $(SP)\ilut utilites $(F)
OBJ = $(addsuffix .obj, $(SRC) )
.SUFFIXES: .obj .mod .f90 .F90 .f
%.obj: %.f90
	$(F90) -o $@ -c $< $(FCOPT1)  
%.obj: %.f
	$(F90) -o $@ -c $< $(FCOPT3) 
$(TARGET): $(OBJ)
	$(F90)  -o $@  $^  $(Lsw)

endif

# Cleaning everything
clean:
ifeq ($(F90),gfortran) 
	del *.o; del *.mod; del *.exe
endif
ifeq ($(F90),ifx) 
	del *.obj; del *.mod; del *.exe
endif
ifeq ($(F90),ifort) 
	del *.obj; del *.mod; del *.exe
endif

cleanout:
	rmdir /s /Q out

