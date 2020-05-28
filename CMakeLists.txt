cmake_minimum_required (VERSION 3.0)
project (HeatModernFortran)
enable_language(Fortran)

# Check whether MPI support is available.
find_package(MPI)
if (MPI_FOUND)
  message(STATUS "modules = ${MPI_Fortran_HAVE_F90_MODULE} ${MPI_Fortran_HAVE_F08_MODULE} ${MPI_Fortran_MODULE_DIR} ")
  message(STATUS "MPI found")
  include_directories(${MPI_INCLUDE_PATH})
  include_directories(${MPI_Fortran_MODULE_DIR})

  # par_error example
  add_library(par_error SHARED par_error.f90)
  target_link_libraries(par_error ${MPI_Fortran_LIBRARIES})
  install(TARGETS par_error DESTINATION ${CMAKE_SOURCE_DIR})

else(MPI_FOUND)
    message(FATAL_ERROR "MPI not found; you NEED it to compile the project")
endif(MPI_FOUND)
