
cmake_minimum_required(VERSION 3.5)

project(3dim CXX)

# Add module path in case this is project root
if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake-modules/" ${CMAKE_MODULE_PATH})
endif()

find_package(Corrade REQUIRED Main)
find_package(Magnum REQUIRED
    DebugTools
    GL
    MeshTools
    Primitives
    Shaders
    Sdl2Application)

find_package(TBB REQUIRED)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_executable(test-particles
    atom.h
    atom.cpp
    constants.h
    logger.h
    logger.cpp
    particle.h
    particle.cpp
    test_particles.cpp
    thread_pool.h
    )

add_executable(3dim WIN32
    3dim.cpp
    atom.h
    atom.cpp
    constants.h
    logger.h
    logger.cpp
    particle.h
    particle.cpp
    arcball/ArcBall.cpp
    thread_pool.h
)

target_link_libraries(3dim PRIVATE
    Corrade::Main
    Magnum::Application
    Magnum::DebugTools
    Magnum::GL
    Magnum::Magnum
    Magnum::MeshTools
    Magnum::Primitives
    Magnum::Shaders
    TBB::tbb)

install(TARGETS 3dim DESTINATION ${MAGNUM_BINARY_INSTALL_DIR})

# Make the executable a default target to build & run in Visual Studio
set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT 3dim)

if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU|Clang)")
  add_compile_options(-g)
  add_compile_options(-std=gnu++20)
endif()
add_compile_options(-std=gnu++20)
