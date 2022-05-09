cmake_minimum_required(VERSION 3.16)

project(AlgorithmEngineering VERSION 1.0 LANGUAGES C CXX)

# Default to development builds.
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Development" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Development" "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif ()

# In development mode,
#  - all warnings are enabled,
#  - full optimization is enabled,
#  - debug information is emitted, and
#  - NDEBUG is NOT defined (asserts are enabled).
if (CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
    set(dev_flags "-Wall -Wextra -Wpedantic -Wnull-dereference -g -O3 -march=native")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(dev_flags "/Wall /permissive- /Zi /O2")
elseif (CMAKE_BUILD_TYPE STREQUAL "Development")
    message(WARNING "Unsupported compiler '${CMAKE_CXX_COMPILER_ID}', development flags NOT set.")
endif ()

set(CMAKE_CXX_FLAGS_DEVELOPMENT "${dev_flags}")
set(CMAKE_C_FLAGS_DEVELOPMENT "${dev_flags}")
set(CMAKE_EXE_LINKER_FLAGS_DEVELOPMENT "${dev_flags}")
set(CMAKE_SHARED_LINKER_FLAGS_DEVELOPMENT "${dev_flags}")
unset(dev_flags)
mark_as_advanced(
    CMAKE_CXX_FLAGS_DEVELOPMENT
    CMAKE_C_FLAGS_DEVELOPMENT
    CMAKE_EXE_LINKER_FLAGS_DEVELOPMENT
    CMAKE_SHARED_LINKER_FLAGS_DEVELOPMENT)

# Enable link time optimization (for all build types).
##include(CheckIPOSupported)
##check_ipo_supported(RESULT result)
##if (result)
##    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
##endif ()
##unset(result)

# Require C++17 standard conformance.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

# Multi-threading support.
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads)

add_library(base_target INTERFACE)
target_include_directories(base_target INTERFACE ${CMAKE_SOURCE_DIR})

add_subdirectory(src)
add_subdirectory(benchmark)

message(${CMAKE_SYSTEM})

