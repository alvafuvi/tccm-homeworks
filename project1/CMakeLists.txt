cmake_minimum_required(VERSION 3.5.0)
project(MP2 VERSION 0.1.0 LANGUAGES C)


set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

set(SOURCE_FILE ${SOURCE_DIR}/mp2.c)

include_directories(/usr/local/include)

link_directories(/usr/local/lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SOURCE_DIR})

add_executable(mp2 ${SOURCE_FILE})

target_link_libraries(mp2 trexio)

