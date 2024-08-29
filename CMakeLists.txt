cmake_minimum_required(VERSION 3.10)
project(simpar)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(Boost 1.66.0 REQUIRED)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CONFIG_DIR "${CMAKE_SOURCE_DIR}/config")
set(OUTPUT_DIR "${CMAKE_BINARY_DIR}/bin")
set(CONFIG_FILES "Big_Network.json"
                "Network.json"
                "Edge_1.json"
                "Edge_2.json"
                "Edge_3.json"
                "Edge_4.json"
                "Edge_5.json"
                "Edge_6.json"
                "Edge_7.json"
                "Edge_8.json"
)

# Iterate over each file and create a symlink
foreach(file ${CONFIG_FILES})
    add_custom_command(
        OUTPUT "${OUTPUT_DIR}/${file}"
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${CONFIG_DIR}/${file}" "${OUTPUT_DIR}/${file}"
        DEPENDS "${CONFIG_DIR}/${file}"
        COMMENT "Creating symlink for ${file}"
    )

    list(APPEND SYMLINKS "${OUTPUT_DIR}/${file}")
endforeach()

add_custom_target(simpar ALL DEPENDS ${SYMLINKS})

# Optionally, make the actual build target depend on the symlink target
add_dependencies(${PROJECT_NAME} simpar)

add_executable(sequential sequential.cpp )
add_executable(network network.cpp network_debug.cpp)
add_executable(network_par network_par.cpp network_par_debug.cpp)
add_executable(network_hybrid network_hybrid.cpp)

option(ENABLE_DUMP "Enable dump" OFF)
set(VIS_DUMP_OPT "DEBUG")

if(ENABLE_DUMP)
    target_compile_definitions(sequential       PUBLIC ${VIS_DUMP_OPT})
    target_compile_definitions(network          PUBLIC ${VIS_DUMP_OPT})
    target_compile_definitions(network_par      PUBLIC ${VIS_DUMP_OPT})
    target_compile_definitions(network_hybrid   PUBLIC ${VIS_DUMP_OPT})
endif()

target_compile_options(sequential PUBLIC 
    -Wall -Wextra -pedantic -Wunused-parameter)

target_compile_options(network PUBLIC 
    -Wall -Wextra -pedantic -Wunused-parameter)

target_compile_options(network_par PUBLIC 
    -Wall -Wextra -pedantic -Wunused-parameter)
    
target_compile_options(network_hybrid PUBLIC 
    -Wall -Wextra -pedantic -Wunused-parameter)

target_link_libraries(network_hybrid PRIVATE OpenMP::OpenMP_CXX)
