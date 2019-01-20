cmake_minimum_required(VERSION 3.0)

project(trw_s)

set(CMAKE_BUILD_TYPE "Debug")

aux_source_directory(. DIR_SRCS_)

add_executable(trw_s ${DIR_SRCS_} typePotts2.h)

install(TARGETS trw_s RUNTIME DESTINATION bin)
