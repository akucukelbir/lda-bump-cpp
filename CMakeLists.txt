cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8.0)

project(lda_hdp_bumping)
set(CMAKE_BUILD_TYPE Release)

add_definitions(-std=c++11)
add_definitions(-stdlib=libc++)

include_directories(/usr/local/include/eigen3/)

find_package(Boost COMPONENTS filesystem system)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(docopt)
include_directories(util)
include_directories(lda)

set(cc_sources docopt.cpp driver.cpp)
add_executable(driver ${cc_sources})
target_link_libraries(driver ${Boost_LIBRARIES})
