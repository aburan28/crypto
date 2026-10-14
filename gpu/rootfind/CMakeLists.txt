# gpu/rootfind -- CMake entry point (same pattern as gpu/mlwe: CUDA is
# optional and detected with check_language).
#
#   cmake -S gpu/rootfind -B build-rootfind && cmake --build build-rootfind && ctest --test-dir build-rootfind
#
# Without nvcc:   librootfind_oracle.a (host emulation) + rootfind_test_cpu
# With nvcc:      additionally librootfind_cuda.a + rootfind_bench, built for
#                 CMAKE_CUDA_ARCHITECTURES (default "120" = RTX PRO 6000
#                 Blackwell; CMake emits -gencode arch=compute_120,code=sm_120
#                 and the PTX fallback for "120-real;120-virtual").
cmake_minimum_required(VERSION 3.24)
project(rootfind_oracle LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Threads REQUIRED)

set(RF_HDRS fp64.cuh fp256x64.cuh rootfind.cuh rootfind_host.hpp rootfind_api.h)

# Emulation-only library: what the Rust build (rs/build.rs) also compiles.
add_library(rootfind_oracle STATIC rootfind_api.cpp ${RF_HDRS})
target_compile_definitions(rootfind_oracle PRIVATE RF_HAVE_CUDA=0)
target_include_directories(rootfind_oracle PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(rootfind_oracle PUBLIC Threads::Threads)

add_executable(rootfind_test_cpu test_cpu.cpp)
target_link_libraries(rootfind_test_cpu PRIVATE rootfind_oracle)
enable_testing()
add_test(NAME rootfind_test_cpu COMMAND rootfind_test_cpu)

include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  if(NOT CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES "120-real;120-virtual")
  endif()
  message(STATUS "rootfind: CUDA ${CMAKE_CUDA_COMPILER_VERSION} found, archs ${CMAKE_CUDA_ARCHITECTURES}")
  add_library(rootfind_cuda STATIC rootfind_api.cpp rootfind_cuda.cu ${RF_HDRS})
  target_compile_definitions(rootfind_cuda PRIVATE RF_HAVE_CUDA=1)
  target_include_directories(rootfind_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(rootfind_cuda PUBLIC Threads::Threads)
  set_target_properties(rootfind_cuda PROPERTIES CUDA_STANDARD 17 CUDA_SEPARABLE_COMPILATION OFF)
  target_compile_options(rootfind_cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>)

  add_executable(rootfind_bench bench.cu)
  target_compile_definitions(rootfind_bench PRIVATE RF_HAVE_CUDA=1)
  target_link_libraries(rootfind_bench PRIVATE rootfind_cuda)
  set_target_properties(rootfind_bench PROPERTIES CUDA_STANDARD 17)
  add_test(NAME rootfind_cuda_selftest COMMAND rootfind_bench selftest 61 1024 4096)
  add_test(NAME rootfind_cuda_p256 COMMAND rootfind_bench p256 256 1024)
else()
  message(STATUS "rootfind: no CUDA compiler; building host emulation only")
endif()
