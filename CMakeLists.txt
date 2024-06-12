cmake_minimum_required(VERSION 3.10)
project(simpar)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(sequential sequential.cpp )
add_executable(network network.cpp network_debug.cpp)
add_executable(network_par network_par.cpp network_par_debug.cpp)
add_executable(network_hybrid network_hybrid.cpp)

target_compile_options(sequential PUBLIC 
    -Wall -Wextra -Wunused-parameter)

target_compile_options(network PUBLIC 
    -Wall -Wextra -Wunused-parameter)

target_compile_options(network_par PUBLIC 
    -Wall -Wextra -Wunused-parameter)

target_compile_options(network_hybrid PUBLIC 
    -Wall -Wextra -Wunused-parameter)

target_link_libraries(network_hybrid PRIVATE OpenMP::OpenMP_CXX)
