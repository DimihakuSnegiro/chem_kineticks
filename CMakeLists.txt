cmake_minimum_required(VERSION 3.16)
project(chemar LANGUAGES CXX)

add_executable(chemar main.cpp)
set_target_properties(chemar PROPERTIES CXX_EXTENSIONS OFF)
target_compile_features(chemar PRIVATE cxx_std_17)
target_include_directories(chemar PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party")

if(MSVC)
    target_compile_options(chemar PRIVATE /W4)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(chemar PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(MINGW)
    target_link_options(chemar PRIVATE -municode -static)
endif()
