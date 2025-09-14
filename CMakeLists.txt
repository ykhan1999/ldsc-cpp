cmake_minimum_required(VERSION 3.20)
project(ldsc_cpp LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(LDSC_ENABLE_OPENMP "Enable OpenMP parallelism" ON)
option(LDSC_ENABLE_IPO "Enable LTO/IPO if supported" ON)

# Sources
add_executable(ldsc
  src/main.cpp
  src/munge.cpp
  src/ph2.cpp
  src/logger.cpp
)

# Include directories
target_include_directories(ldsc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Warnings
if (MSVC)
  target_compile_options(ldsc PRIVATE /W4 /permissive-)
else()
  target_compile_options(ldsc PRIVATE -Wall -Wextra -Wpedantic)
endif()

# OpenMP (optional)
if (LDSC_ENABLE_OPENMP)
  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
    target_link_libraries(ldsc PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(ldsc PRIVATE LDSC_HAS_OPENMP=1)
  endif()
endif()

# IPO/LTO
include(CheckIPOSupported)
if (LDSC_ENABLE_IPO)
  check_ipo_supported(RESULT ipo_ok OUTPUT ipo_msg)
  if (ipo_ok)
    set_property(TARGET ldsc PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
endif()

# Install
include(GNUInstallDirs)
install(TARGETS ldsc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Packaging (archive)
set(CPACK_PACKAGE_NAME "ldsc-cpp")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
include(CPack)

