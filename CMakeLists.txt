cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
  cmake_policy(SET CMP0071 NEW)
endif()

# Set the project
PROJECT(TIGER VERSION 1.0 DESCRIPTION "Topological Interlocking GEneratoR" LANGUAGES CXX)

add_subdirectory(src)

