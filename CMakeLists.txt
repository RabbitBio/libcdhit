cmake_minimum_required(VERSION 3.10)
project(libcdhit)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/libcdhit)

set(LIB_NAME libcdhit)

#option(BUILD_SHARED_LIBS "Build shared library" OFF)

file(GLOB_RECURSE SOURCES "${SOURCE_DIR}/*.cpp")

#add_library(${LIB_NAME} STATIC ${SOURCES})
add_library(${LIB_NAME} SHARED ${SOURCES})

if(EXISTS ${INCLUDE_DIR})
    target_include_directories(${LIB_NAME} PUBLIC ${INCLUDE_DIR})
endif()

target_compile_options(${LIB_NAME} PRIVATE -Wall -O3)

install(TARGETS ${LIB_NAME} DESTINATION lib)
install(DIRECTORY ${INCLUDE_DIR}/ DESTINATION include/libcdhit)


