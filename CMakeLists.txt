### CMake project file for myprog
cmake_minimum_required (VERSION 2.6)
project (myprog)
enable_language (Fortran)





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
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -Wall -fimplicit-none -fbounds-check -O -Wuninitialized -ffpe-trap=invalid, zero, overflow -fbacktrace")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g -warn all -check bounds -check uninit -fpe0 -traceback -heap-arrays")
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
