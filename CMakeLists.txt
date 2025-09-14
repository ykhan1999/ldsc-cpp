cmake_minimum_required(VERSION 3.20)
project(ldsc_cpp LANGUAGES CXX VERSION 0.3.0)

# ---- Language/flags ---------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Prefer Release if unspecified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Options
option(LDSC_ENABLE_OPENMP "Enable OpenMP parallelism" ON)
option(LDSC_ENABLE_IPO     "Enable LTO/IPO if supported" ON)
option(LDSC_PORTABLE_LINUX "Link libstdc++/libgcc statically on Linux (more portable)" ON)
option(LDSC_MSVC_STATIC_CRT "Static MSVC runtime on Windows (/MT)" OFF)

# macOS portability: set a modern, broadly compatible target
if(APPLE)
  # Big Sur is the first common baseline across Intel/Apple Silicon
  set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "macOS deployment target")
  # Ensure runtime search paths work well for a standalone binary
  set(CMAKE_MACOSX_RPATH ON)
  set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
  set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
endif()

# Windows portability: allow switching to static CRT easily
if(MSVC)
  cmake_policy(SET CMP0091 NEW) # respect CMAKE_MSVC_RUNTIME_LIBRARY
  if(LDSC_MSVC_STATIC_CRT)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  else()
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL$<$<CONFIG:Debug>:Debug>")
  endif()
endif()

# Linux portability: prefer static libstdc++/libgcc for wider distro coverage
if(UNIX AND NOT APPLE AND LDSC_PORTABLE_LINUX)
  add_link_options(-static-libstdc++ -static-libgcc)
endif()

# Warnings
if (MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# ---- Sources/targets --------------------------------------------------------
add_executable(ldsc
  src/main.cpp
  src/munge.cpp
  src/ph2.cpp
  src/logger.cpp
)

target_compile_definitions(ldsc PRIVATE
  PROJECT_NAME_STR="ldsc-cpp"
  PROJECT_VERSION_STR="${PROJECT_VERSION}"
  PROJECT_LICENSE_STR="GPL-3.0-or-later"
  PROJECT_URL_STR="https://github.com/ykhan1999/ldsc-cpp"
  PROJECT_AUTHOR_STR="Yousef Khan"
)

target_include_directories(ldsc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Project metadata exposed to the binary (for --version / banner)
target_compile_definitions(ldsc PRIVATE
  PROJECT_NAME_STR="ldsc-cpp"
  PROJECT_VERSION_STR="${PROJECT_VERSION}"
  PROJECT_LICENSE_STR="GPL-3.0-or-later"
  PROJECT_URL_STR="https://github.com/ykhan1999/ldsc-cpp"
  PROJECT_AUTHOR_STR="Yousef Khan"
)

# OpenMP (optional; resilient on macOS/Homebrew)
if (LDSC_ENABLE_OPENMP)
  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
    target_link_libraries(ldsc PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(ldsc PRIVATE LDSC_HAS_OPENMP=1)
  elseif(APPLE)
    find_library(OMP_LIB omp HINTS /opt/homebrew/lib /usr/local/lib)
    if(OMP_LIB)
      target_link_libraries(ldsc PRIVATE ${OMP_LIB})
      target_compile_options(ldsc PRIVATE -Xpreprocessor -fopenmp)
      target_compile_definitions(ldsc PRIVATE LDSC_HAS_OPENMP=1)
    endif()
  endif()
endif()

# IPO/LTO (Release/RelWithDebInfo)
include(CheckIPOSupported)
if (LDSC_ENABLE_IPO)
  check_ipo_supported(RESULT ipo_ok OUTPUT ipo_msg)
  if (ipo_ok)
    set_property(TARGET ldsc PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    set_property(TARGET ldsc PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
  else()
    message(STATUS "IPO not supported: ${ipo_msg}")
  endif()
endif()

# Install
include(GNUInstallDirs)
install(TARGETS ldsc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Packaging (CPack)
set(CPACK_PACKAGE_NAME "ldsc-cpp")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_GENERATOR "TGZ;ZIP")
include(CPack)

# ---- Nice build summary -----------------------------------------------------
message(STATUS "======== ldsc_cpp build config ========")
message(STATUS "CMAKE_BUILD_TYPE      : ${CMAKE_BUILD_TYPE}")
message(STATUS "OpenMP enabled        : ${LDSC_ENABLE_OPENMP}")
message(STATUS "IPO/LTO enabled       : ${LDSC_ENABLE_IPO}")
if(APPLE)
  message(STATUS "macOS target          : ${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif()
if(WIN32)
  message(STATUS "MSVC static CRT       : ${LDSC_MSVC_STATIC_CRT}")
endif()
if(UNIX AND NOT APPLE)
  message(STATUS "Linux static libstdc++: ${LDSC_PORTABLE_LINUX}")
endif()
message(STATUS "=======================================")

