cmake_minimum_required(VERSION 3.3)
project(valiant_parsing)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    README.md main.cpp grammar.cpp grammar.h matrix.h)

add_executable(valiant_parsing ${SOURCE_FILES})