cmake_minimum_required(VERSION 3.5)

project(fixed_point LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(fixed_point main.cpp
    decimal.h
    longtypes.h)

if(MSVC)
  target_compile_options(fixed_point PRIVATE /W4 /WX)
else()
  target_compile_options(fixed_point PRIVATE -Wall -Wextra -Wpedantic -Werror -std=c++17)
endif()


include(GNUInstallDirs)
install(TARGETS fixed_point
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
