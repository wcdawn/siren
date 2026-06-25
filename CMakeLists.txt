cmake_minimum_required(VERSION 3.14)
project(siren VERSION 0.2.0 DESCRIPTION "A one-dimensional PN neutron transport solver." LANGUAGES Fortran)

add_executable(${PROJECT_NAME}.x)

target_sources(${PROJECT_NAME}.x
  PRIVATE
    src/analytic.f90
    src/constant.f90
    src/diffusion_block.f90
    src/diffusion.f90
    src/exception_handler.f90
    src/fileio.f90
    src/geometry.f90
    src/input.f90
    src/kind.f90
    src/linalg.f90
    src/main.f90
    src/numeric.f90
    src/output.f90
    src/power.f90
    src/timer.f90
    src/transport_block.f90
    src/transport.f90
    src/xs.f90
)

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME}.x OpenMP::OpenMP_Fortran)

find_package(LAPACK REQUIRED)
target_link_libraries(${PROJECT_NAME}.x LAPACK::LAPACK)

# Note that the CMake language does not allow setting the Fortran_STANDARD from the build instructions
# and it must be specified via the command-line argument.
# TODO debug flags
target_compile_options(${PROJECT_NAME}.x PRIVATE -Wall -Wextra -O3 -std=f2018 -pedantic)
