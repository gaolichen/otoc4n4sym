cmake_minimum_required (VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

SET(MYDEV D:/dev)
SET(BOOST_ROOT ${MYDEV}/boost_1_70_0)
find_package(Boost 1.70)

project (otoc4n4sym)
SET (CMAKE_BUILD_TYPE Release)
# to use gdb enable the following two flags.
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)

add_subdirectory (src)
