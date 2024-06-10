cmake_minimum_required(VERSION 3.10)
project(simpar)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(sequential.out
    sequential.cpp
)

add_executable(network.out
    network.cpp
    network_debug.cpp
)

add_executable(network_par.out
    network_par.cpp
    network_par_debug.cpp    
)

add_executable(network_hybrid.out
    network_hybrid.cpp    
)
