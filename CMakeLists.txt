cmake_minimum_required(VERSION 3.10)
project(libcdhit)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/libcdhit)

set(LIB_NAME cdhitlib)

option(AVX2 "Enable AVX2 optimizations if available" OFF)
option(AVX512 "Enable AVX512 optimizations if available" OFF)
if(AVX2)
    add_compile_options(-mavx2)
endif()
if(AVX512)
    add_compile_options(-mavx512f -mavx512bw -mavx512vl)
endif()

# ---- Parasail ----
set(PARASAIL_ROOT "" CACHE PATH "Parasail install prefix")

if (PARASAIL_ROOT)
  target_include_directories(${LIB_NAME} PUBLIC "${PARASAIL_ROOT}/include")

  find_library(PARASAIL_LIBRARY
    NAMES parasail
    HINTS "${PARASAIL_ROOT}"
    PATH_SUFFIXES lib lib64
    REQUIRED
  )

  target_link_libraries(${LIB_NAME} PUBLIC "${PARASAIL_LIBRARY}")
endif()
#-------------------

# ---- edlib ----
# 查找 edlib 动态库（通过 -DEDLIB_ROOT 指定或自动查找）
if(EDLIB_ROOT)
    set(EDLIB_INCLUDE_DIR ${EDLIB_ROOT}/include)
    find_library(EDLIB_LIBRARY NAMES edlib libedlib.so PATHS ${EDLIB_ROOT} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
    find_path(EDLIB_INCLUDE_DIR edlib.h PATHS /usr/include /usr/local/include)
    find_library(EDLIB_LIBRARY NAMES edlib libedlib.so PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)
endif()

if(NOT EDLIB_INCLUDE_DIR OR NOT EDLIB_LIBRARY)
    message(FATAL_ERROR "edlib not found. Please specify -DEDLIB_ROOT=<path> or install edlib system-wide.")
endif()
#-------------------

#option(BUILD_SHARED_LIBS "Build shared library" OFF)

file(GLOB_RECURSE SOURCES "${SOURCE_DIR}/*.cpp")

#add_library(${LIB_NAME} STATIC ${SOURCES})
add_library(${LIB_NAME} SHARED ${SOURCES})
target_link_libraries(${LIB_NAME} PUBLIC OpenMP::OpenMP_CXX)  # Add this line to link OpenMP to the library
target_include_directories(${LIB_NAME} PUBLIC ${EDLIB_INCLUDE_DIR})
target_link_libraries(${LIB_NAME} PRIVATE ${EDLIB_LIBRARY})


if(EXISTS ${INCLUDE_DIR})
    target_include_directories(${LIB_NAME} PUBLIC ${INCLUDE_DIR})
endif()

target_compile_options(${LIB_NAME} PRIVATE -Wall -O3)

add_executable(test_cdhit tests/test.cpp)
target_include_directories(test_cdhit PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_cdhit PRIVATE cdhitlib OpenMP::OpenMP_CXX)

add_executable(test_cc tests/test_cc.cpp)
target_include_directories(test_cc PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_cc PRIVATE z cdhitlib OpenMP::OpenMP_CXX)

install(TARGETS ${LIB_NAME} DESTINATION lib)
install(DIRECTORY ${INCLUDE_DIR}/ DESTINATION include/libcdhit)


