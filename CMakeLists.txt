<<<<<<< HEAD
### CMake project file for myprog 
cmake_minimum_required (VERSION 3.1.1)
=======
### CMake project file for myprog
cmake_minimum_required (VERSION 2.6)
<<<<<<< HEAD
>>>>>>> 9f36052... new code
=======
>>>>>>> 9f36052... new code
project (myprog)
enable_language (Fortran) 





<<<<<<< HEAD
<<<<<<< HEAD
add_executable (a.out main.f90 global.f90 hamiltonian.f90 basis.f90 boundary.f90 inner.f90 outer.f90 asymptote.f90) 

### fgsl library 
include_directories (/usr/local/include/fgsl)
target_link_libraries (a.out -lfgsl)

### lapack library 
# target_link_libraries (a.out -llapack) 
include_directories (-i8)
include_directories (/opt/intel/mkl/include)
include_directories (/opt/intel/mkl/include/intel64/ilp64) 
target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_blas95_ilp64.a)  
target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_lapack95_ilp64.a)
target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_intel_ilp64.a)
target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_core.a) 
target_link_libraries (a.out /opt/intel/mkl/lib/libmkl_sequential.a) 
target_link_libraries (a.out -lpthread)
target_link_libraries (a.out -lm) 

### my library 
include_directories (/Users/hyuna917/Dropbox/Library/Fortran/include/) 
target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/libconst.a) 
# target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/libfft.a) 
target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/liblinear.a) 
# target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/libplot.a) 
# target_link_libraries (a.out /Users/hyuna917/Dropbox/Library/Fortran/lib/libsde.a) 





### make sure that the default is a RELEASE
=======
=======
>>>>>>> 9f36052... new code
add_executable (a.out
    main.f90
    global.f90
    hamiltonian.f90
    inner.f90
    outer.f90
    )
target_link_libraries (a.out
    fgsl
    # lapack
    myconst.150401
    # myfft
    mylinear.150326
    # myplot
    # mysde
    )

### directories
include_directories (/usr/include /usr/local/include $ENV{HOME}/include)
link_directories (/usr/lib /usr/local/lib $ENV{HOME}/lib)
target_link_libraries (a.out -I$ENV{HOME}/include -L$ENV{HOME}/lib)

### intel library
# if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#     set(MKL_INCLUDE /opt/intel/mkl/include)
#     set(MKL_LIBRARY /opt/intel/mkl/lib)
# elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#     set(MKL_INCLUDE /opt/intel/mkl/include)
#     set(MKL_LIBRARY /opt/intel/mkl/lib/intel64)
# endif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# include_directories (${MKL_INCLUDE})
# link_directories (${MKL_LIBRARY})
# target_link_libraries (a.out ${MKL_LIBRARY}/libmkl_blas95_ilp64.a)
# target_link_libraries (a.out ${MKL_LIBRARY}/libmkl_lapack95_ilp64.a)
# target_link_libraries (a.out ${MKL_LIBRARY}/libmkl_intel_ilp64.a)
# target_link_libraries (a.out ${MKL_LIBRARY}/libmkl_core.a)
# target_link_libraries (a.out ${MKL_LIBRARY}/libmkl_sequential.a)
# target_link_libraries (a.out -lpthread)
# target_link_libraries (a.out -lm)





# make sure that the default is a RELEASE
<<<<<<< HEAD
>>>>>>> 9f36052... new code
=======
>>>>>>> 9f36052... new code
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the NORMA directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
<<<<<<< HEAD
<<<<<<< HEAD
	### gfortran
	message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
	set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
	### ifort (untested)
	message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
	set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
=======
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -Wall -fimplicit-none -fbounds-check -O -Wuninitialized -ffpe-trap=invalid, zero, overflow -fbacktrace")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
=======
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -Wall -fimplicit-none -fbounds-check -O -Wuninitialized -ffpe-trap=invalid, zero, overflow -fbacktrace")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
>>>>>>> 9f36052... new code
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g -warn all -check bounds -check uninit -fpe0 -traceback -heap-arrays")
<<<<<<< HEAD
>>>>>>> 9f36052... new code
=======
>>>>>>> 9f36052... new code
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
