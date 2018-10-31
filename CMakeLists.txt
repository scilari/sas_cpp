cmake_minimum_required(VERSION 3.12)
project(sas_cpp)

set(CMAKE_CXX_STANDARD 11)

add_executable(sas_cpp main.cpp)

include_directories(${MY_SOURCE_DIR}/src)

link_directories(${MY_SOURCE_DIR}/src)

target_link_libraries(sas_cpp sas)

add_subdirectory(src)

add_subdirectory(test)

enable_testing()

