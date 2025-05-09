cmake_minimum_required(VERSION 3.16)
project(discharge CXX)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Boost REQUIRED COMPONENTS thread)
find_package(spdlog REQUIRED)

add_executable(a.out main.cpp near_triangulation.cpp cartwheel.cpp configuration.cpp rule.cpp basewheel.cpp)
target_compile_options(a.out PUBLIC -O2 -Wall -Wextra)
target_compile_features(a.out PUBLIC cxx_std_20)
target_link_libraries(a.out PRIVATE 
    Boost::boost Boost::program_options
    spdlog::spdlog)

add_executable(send send.cpp near_triangulation.cpp cartwheel.cpp configuration.cpp rule.cpp basewheel.cpp)
target_compile_options(send PUBLIC -O2 -Wall -Wextra)
target_compile_features(send PUBLIC cxx_std_20)
target_link_libraries(send PRIVATE 
    Boost::boost Boost::program_options
    spdlog::spdlog)

add_executable(unique unique.cpp near_triangulation.cpp cartwheel.cpp configuration.cpp rule.cpp basewheel.cpp)
target_compile_options(unique PUBLIC -O2 -Wall -Wextra)
target_compile_features(unique PUBLIC cxx_std_20)
target_link_libraries(unique PRIVATE 
    Boost::boost Boost::program_options
    spdlog::spdlog)

add_executable(tile tile.cpp near_triangulation.cpp cartwheel.cpp configuration.cpp rule.cpp basewheel.cpp)
target_compile_options(tile PUBLIC -O2 -Wall -Wextra)
target_compile_features(tile PUBLIC cxx_std_20)
target_link_libraries(tile PRIVATE 
    Boost::boost Boost::program_options Boost::thread
    spdlog::spdlog)
