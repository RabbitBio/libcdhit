cmake_minimum_required(VERSION 3.10)
project(libcdhit LANGUAGES CXX)

# ===== OpenMP =====
find_package(OpenMP REQUIRED)

# ===== C++ 标准 =====
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 默认 Release（若未显式指定）
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# 可选：是否启用 SIMD（AVX2/POPCNT/BMI）优化
option(ENABLE_SIMD "Enable AVX2/POPCNT/BMI compiler flags when available" ON)
# 可选：Link Time Optimization（IPO/LTO）
option(ENABLE_LTO "Enable Interprocedural Optimization (LTO) if supported" ON)

set(SOURCE_DIR  ${CMAKE_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/libcdhit)
set(LIB_NAME    cdhitlib)

file(GLOB_RECURSE SOURCES "${SOURCE_DIR}/*.cpp")

# ===== 库目标 =====
add_library(${LIB_NAME} SHARED ${SOURCES})
target_link_libraries(${LIB_NAME} PUBLIC OpenMP::OpenMP_CXX)

if(EXISTS ${INCLUDE_DIR})
    target_include_directories(${LIB_NAME} PUBLIC ${INCLUDE_DIR})
endif()

# 基础警告 & 优化
if (MSVC)
  target_compile_options(${LIB_NAME} PRIVATE /O2 /W3)
else()
  target_compile_options(${LIB_NAME} PRIVATE -O3 -Wall -g -Wextra)
endif()

# ===== SIMD（AVX2/POPCNT/BMI）按编译器可用性追加 =====
include(CheckCXXCompilerFlag)

if (ENABLE_SIMD)
  if (MSVC)
    # MSVC：/arch:AVX2 即可；POPCNT/BMI 通常随 /arch:AVX2 可用
    target_compile_options(${LIB_NAME} PRIVATE /arch:AVX2)
    target_compile_definitions(${LIB_NAME} PRIVATE HAVE_AVX2=1)
  else()
    set(_simd_flags -mavx2 -mpopcnt -mbmi)
    # set(_simd_flags -mavx512f -mavx512bw -mavx512vl -mavx512dq -mavx512cd -mpopcnt -mbmi)
    foreach(f IN LISTS _simd_flags)
      string(REPLACE "-" "_" fvar ${f})
      check_cxx_compiler_flag(${f} HAVE_${fvar})
      if(HAVE_${fvar})
        target_compile_options(${LIB_NAME} PRIVATE ${f})
        if(f STREQUAL "-mavx2")
          target_compile_definitions(${LIB_NAME} PRIVATE HAVE_AVX2=1)
        elseif(f STREQUAL "-mbmi")
          target_compile_definitions(${LIB_NAME} PRIVATE HAVE_BMI=1)
        endif()
      endif()
    endforeach()
  endif()
endif()

# ===== LTO/IPO =====
if (ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_msg)
  if(_ipo_ok)
    set_property(TARGET ${LIB_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
endif()

# ===== 测试可执行文件 =====
add_executable(test_cdhit tests/test.cpp)
target_include_directories(test_cdhit PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_cdhit PRIVATE ${LIB_NAME} OpenMP::OpenMP_CXX)

add_executable(test_cc tests/test_cc.cpp)
target_include_directories(test_cc PRIVATE ${CMAKE_SOURCE_DIR}/include)
# 你原来直接链接 z；这里保持不变（若需更可移植可改为 find_package(ZLIB)）
target_link_libraries(test_cc PRIVATE z ${LIB_NAME} OpenMP::OpenMP_CXX)

# 将相同的优化/指令集选项应用到两个测试目标
if (MSVC)
  target_compile_options(test_cdhit PRIVATE /O2 /W3)
  target_compile_options(test_cc    PRIVATE /O2 /W3)
  if (ENABLE_SIMD)
    target_compile_options(test_cdhit PRIVATE /arch:AVX2)
    target_compile_options(test_cc    PRIVATE /arch:AVX2)
    target_compile_definitions(test_cdhit PRIVATE HAVE_AVX2=1)
    target_compile_definitions(test_cc    PRIVATE HAVE_AVX2=1)
  endif()
else()
  target_compile_options(test_cdhit PRIVATE -O3 -Wall -Wextra)
  target_compile_options(test_cc    PRIVATE -O3 -Wall -Wextra)
  if (ENABLE_SIMD)
    foreach(f IN LISTS _simd_flags)
      string(REPLACE "-" "_" fvar ${f})
      if(DEFINED HAVE_${fvar} AND HAVE_${fvar})
        target_compile_options(test_cdhit PRIVATE ${f})
        target_compile_options(test_cc    PRIVATE ${f})
        if(f STREQUAL "-mavx2")
          target_compile_definitions(test_cdhit PRIVATE HAVE_AVX2=1)
          target_compile_definitions(test_cc    PRIVATE HAVE_AVX2=1)
        elseif(f STREQUAL "-mbmi")
          target_compile_definitions(test_cdhit PRIVATE HAVE_BMI=1)
          target_compile_definitions(test_cc    PRIVATE HAVE_BMI=1)
        endif()
      endif()
    endforeach()
  endif()
endif()

# LTO 也给测试程序
if (ENABLE_LTO AND _ipo_ok)
  set_property(TARGET test_cdhit PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET test_cc    PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# ===== 安装 =====
install(TARGETS ${LIB_NAME} DESTINATION lib)
install(DIRECTORY ${INCLUDE_DIR}/ DESTINATION include/libcdhit)
