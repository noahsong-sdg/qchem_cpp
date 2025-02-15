cmake_minimum_required(VERSION 3.10)
project(HartreeFock)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add Eigen
include_directories(${PROJECT_SOURCE_DIR}/Eigen)
include_directories(${PROJECT_SOURCE_DIR}/include)

# Create library
add_library(hf_lib 
    src/hf_fns.cpp
    src/hartree_fock.cpp
)

# Create executable
add_executable(hf_program src/app.cpp)
target_link_libraries(hf_program PRIVATE hf_lib)
