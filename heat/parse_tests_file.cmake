#=============================================================
# parse_tests_file
#=============================================================

macro(parse_tests_file filename libs nb_procs)
    if (${nb_procs} GREATER 1)
        set(use_mpi TRUE)
    endif()
    get_filename_component(filename_we  ${filename} NAME_WLE)
    set(file_wrap ${CMAKE_CURRENT_BINARY_DIR}/${filename_we}_wrap.F90)
    set(test_names "")
    # write the header of the wrapper file
    file(WRITE  ${file_wrap}
"program ${filename_we}_wrap 
    use iso_fortran_env
    use ${filename_we}
    character(len=256) :: first_arg

    call get_command_argument(1, first_arg)

    select case( trim(first_arg))
")

    
    # now parse the sources to add the tests from each file
    file(READ "${filename}" contents)
    string(REPLACE "\n" ";" rows ${contents})
    foreach(line ${rows})
      string(REGEX MATCHALL "^ *subroutine *test_([A-Za-z_0-9])+ *" found_suite ${line})
      string(REGEX REPLACE "^ *subroutine *test_(([A-Za-z_0-9])+) *" "\\1" suite_name "${found_suite}")
      foreach(i ${suite_name})
      file(APPEND ${file_wrap}
"        case('${suite_name}')
              call test_${suite_name}
")
      set(test_names ${test_names} ${suite_name})
      endforeach()
    endforeach()
    # write the end of the wrapper
    file(APPEND ${file_wrap}
"        case default
     end select
end program ${filename_we}_wrap")
add_executable(${filename_we} ${filename_we}_wrap.F90  ${filename} asserts.f90)
if (use_mpi)
    target_include_directories(${filename_we} PRIVATE ${MPI_INCLUDE_PATH})
    #target_compile_options(${filename_we} PRIVATE ${Coarray_COMPILE_OPTIONS})
    target_link_libraries(${filename_we} PRIVATE "${libs}" ${MPI_Fortran_LIBRARIES})
else()
    target_link_libraries(${filename_we} PRIVATE "${libs}")
endif()
enable_testing()
foreach(name ${test_names})
    if (use_mpi)
        add_test(
        NAME ${name} 
        COMMAND ${MPIEXEC} ${ADD_OVERSUBSCRIBE_TO_MPI} ${MPIEXEC_NUMPROC_FLAG} ${nb_procs} ${filename_we} ${name}
        ) 
    else()
        add_test(NAME ${name} COMMAND ${filename_we} ${name})
    endif()
endforeach()
endmacro(parse_tests_file)
