cmake_minimum_required(VERSION 3.16)
project(PlasticMuonTiming)

set(CMAKE_CXX_STANDARD 17)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(plasticMuon main.cc ${SOURCES} ${HEADERS})
target_link_libraries(plasticMuon ${Geant4_LIBRARIES})

set(RUN_DIR ${PROJECT_SOURCE_DIR}/run)
file(MAKE_DIRECTORY ${RUN_DIR})

# Copy executable after compilation
add_custom_command(TARGET plasticMuon POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            $<TARGET_FILE:plasticMuon>
            ${RUN_DIR}/plasticMuon
)

execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            ../macros
            ${RUN_DIR}/macros
)

execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            ../analysis
            ${RUN_DIR}/analysis
)
