cmake_minimum_required(VERSION 3.16)
project(ips4o_cli_bench LANGUAGES CXX)

# --- 基本编译设置 ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release CACHE STRING "")

# --- 选项（可按需 -D 开关） ---
option(OPTIMIZE_NATIVE "为本机CPU优化（-O3 -march=native）" ON)
option(ENABLE_PARALLEL "启用并行路径（代码里会定义 USE_IPS4O_PARALLEL）" ON)
option(ENABLE_OPENMP   "若可用则链接 OpenMP" OFF)
option(LINK_ATOMIC     "如遇 __atomic_compare_exchange_16 链接错误则开" OFF)

# --- 让 ips4o 当子工程编 ---
# 注意：在 add_subdirectory 之前设置 ips4o 的 CMake 变量更生效
# 如果你不想本机优化，可把下面两个设成 OFF
set(IPS4O_OPTIMIZE_FOR_NATIVE ${OPTIMIZE_NATIVE} CACHE BOOL "" FORCE)
set(IPS4O_USE_MCX16           ${OPTIMIZE_NATIVE} CACHE BOOL "" FORCE)

add_subdirectory(thirdparty/ips4o ips4o-build)

# --- 你的可执行程序 ---
# 改成你的源文件实际位置：src/main.cpp 或 main.cpp
add_executable(bench src/main.cpp)

# 把你自己的 include/ 加进来，这里有 CLI.hpp
target_include_directories(bench PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# 链接 ips4o 目标（由上面的 add_subdirectory 提供）
target_link_libraries(bench PRIVATE ips4o)

# --- 并行与依赖（可选） ---
if(ENABLE_PARALLEL)
  target_compile_definitions(bench PRIVATE USE_IPS4O_PARALLEL)
  # 如果系统已安装 TBB，这里尝试链接；找不到也没关系，ips4o 子项目会决定是否需要它
  find_package(TBB QUIET)
  if(TBB_FOUND)
    target_link_libraries(bench PRIVATE TBB::tbb)
  endif()
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(bench PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()

# --- 编译优化（不想改全局 CFLAGS 的就只对目标加） ---
if(OPTIMIZE_NATIVE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bench PRIVATE -O3 -DNDEBUG -march=native)
  elseif(MSVC)
    target_compile_options(bench PRIVATE /O2)
  endif()
endif()

# --- 部分平台需要手动链接 libatomic（只在报错时开） ---
if(LINK_ATOMIC)
  if(UNIX AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_link_libraries(bench PRIVATE atomic)
  endif()
endif()
