cmake_minimum_required(VERSION 3.7)
project(hydra_pp CXX)

set(CMAKE_CXX_STANDARD 14)

set(POST_PROC
        HydraPostProcessAPI.h
        HydraPostProcessAPI.cpp
        HydraPostProcessCommon.h
        HydraPostProcessCommon.cpp
        HydraPostProcessSpecial.h
        HydraPostProcessSpecial.cpp
        NonLocalMeans.cpp)

set(HYDRA_API_SRC_DIR ${HYDRA_API_ROOT_DIR}/hydra_api)

MESSAGE(WARNING "${HYDRA_API_SRC_DIR}")

include_directories(${HYDRA_API_SRC_DIR})

add_library(hydra_pp STATIC ${POST_PROC})



