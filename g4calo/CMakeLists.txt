cmake_minimum_required(VERSION 3.16)
project(g4calo CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 11 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
add_executable(calo calo.cc ${sources})
target_include_directories(calo PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(calo ${Geant4_LIBRARIES})
