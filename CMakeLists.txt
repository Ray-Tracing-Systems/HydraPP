cmake_minimum_required(VERSION 3.7)
project(hydra_pp CXX)

add_subdirectory (tests)

set(CMAKE_CXX_STANDARD 14)

set(POST_PROC
        HydraPostProcessAPI.h
        HydraPostProcessAPI.cpp
        HydraPostProcessCommon.h
        HydraPostProcessCommon.cpp
        HydraPostProcessSpecial.h
        HydraPostProcessSpecial.cpp
        NonLocalMeans.cpp)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
                
set(HYDRA_API_ROOT_DIR ${CMAKE_SOURCE_DIR}/../HydraAPI)
set(HYDRA_API_SRC_DIR ${HYDRA_API_ROOT_DIR}/hydra_api)
set(HYDRA_API_LIB_DIR ${HYDRA_API_ROOT_DIR}/bin)
set(HYDRA_API_UTILS_SRC_DIR ${HYDRA_API_ROOT_DIR}/utils)

#MESSAGE(WARNING "${HYDRA_API_ROOT_DIR}")

include_directories(${HYDRA_API_SRC_DIR})

add_library(hydra_pp STATIC ${POST_PROC})



