cmake_minimum_required(VERSION 3.1.3)
project(Mecacell)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O0 -std=c++14 -g -Wall -Wextra -Wundef -Wshadow -Winit-self -Wconversion -Wformat=2 -pedantic ")
add_subdirectory(mecacell)
add_subdirectory(mecacellviewer)
add_subdirectory(tests)
