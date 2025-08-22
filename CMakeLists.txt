cmake_minimum_required(VERSION 3.10)
project(libcdhit)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/libcdhit)

set(LIB_NAME cdhitlib)

#option(BUILD_SHARED_LIBS "Build shared library" OFF)

file(GLOB_RECURSE SOURCES "${SOURCE_DIR}/*.cpp")

#add_library(${LIB_NAME} STATIC ${SOURCES})
add_library(${LIB_NAME} SHARED ${SOURCES})
target_link_libraries(${LIB_NAME} PUBLIC OpenMP::OpenMP_CXX)  # Add this line to link OpenMP to the library

if(EXISTS ${INCLUDE_DIR})
    target_include_directories(${LIB_NAME} PUBLIC ${INCLUDE_DIR})
endif()

target_compile_options(${LIB_NAME} PRIVATE -Wall -O3)

add_executable(test_cdhit tests/test.cpp)
target_include_directories(test_cdhit PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_cdhit PRIVATE cdhitlib OpenMP::OpenMP_CXX)

add_executable(test_cc tests/test_cc.cpp)
target_include_directories(test_cc PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_cc PRIVATE cdhitlib OpenMP::OpenMP_CXX)

install(TARGETS ${LIB_NAME} DESTINATION lib)
install(DIRECTORY ${INCLUDE_DIR}/ DESTINATION include/libcdhit)


