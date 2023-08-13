cmake_minimum_required(VERSION 3.5)

enable_language(Fortran)

project(droplet Fortran)

option(use_OpenMP "enable OpenMP directives" OFF)

if(use_OpenMP)
    find_package(OpenMP REQUIRED)
endif()

#コンパイルオプション指定用関数
function(SETCOMPILEOPTIONS name)

    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        target_compile_options(${name} PRIVATE
            $<$<CONFIG:Debug>:-O0 -traceback > # -warn all -check all 
            $<$<CONFIG:Release>:-O2 -fPIC>
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        target_compile_options(${name} PRIVATE
            $<$<CONFIG:Debug>: -O0 -g -fcheck=bounds>
            $<$<CONFIG:Release>: -O2>
        )
    endif()

    if(OpenMP_Fortran_FOUND)
        target_compile_options(${name} PRIVATE 
            "${OpenMP_Fortran_FLAGS}"
        )
        target_link_libraries(${name}
            OpenMP::OpenMP_Fortran 
        )
    endif()

endfunction(SETCOMPILEOPTIONS)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

add_subdirectory(lib)
add_subdirectory(src)
