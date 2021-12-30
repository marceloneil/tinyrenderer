cmake_minimum_required (VERSION 3.10)
project (tinyrenderer)

include(CheckCXXCompilerFlag)

function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()

enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Wextra")
enable_cxx_compiler_flag_if_supported("-pedantic")
enable_cxx_compiler_flag_if_supported("-std=c++14")
enable_cxx_compiler_flag_if_supported("-O3")
enable_cxx_compiler_flag_if_supported("-fopenmp")

# Method used to rasterize triangles
if(NOT TRIANGLE_METHOD)
    set(TRIANGLE_METHOD BARYCENTRIC)
endif()

# Method used to compute barycentric coordinates
if(NOT BARYCENTRIC_METHOD)
    set(BARYCENTRIC_METHOD SIMPLE)
endif()

add_definitions(-DTRIANGLE_${TRIANGLE_METHOD})
add_definitions(-DBARYCENTRIC_${BARYCENTRIC_METHOD})

file(GLOB SOURCES *.h *.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

