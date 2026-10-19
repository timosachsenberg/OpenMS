cmake_minimum_required(VERSION 3.24)
project(OpenMSXMLParserBenchmark LANGUAGES CXX)

find_package(OpenMS CONFIG REQUIRED)
find_package(Threads REQUIRED)
add_executable(xml_parser_benchmark benchmark.cpp)
target_compile_features(xml_parser_benchmark PRIVATE cxx_std_23)
target_link_libraries(xml_parser_benchmark PRIVATE OpenMS::OpenMS Threads::Threads ${CMAKE_DL_LIBS})
add_executable(xml_schema_validation_probe schema_validation_probe.cpp)
target_compile_features(xml_schema_validation_probe PRIVATE cxx_std_23)
target_link_libraries(xml_schema_validation_probe PRIVATE OpenMS::OpenMS)
if(UNIX)
  add_executable(xml_benchmark_launcher launcher.cpp)
  target_compile_features(xml_benchmark_launcher PRIVATE cxx_std_17)
endif()

enable_testing()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME xml_benchmark_smoke
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_benchmark.py
    $<TARGET_FILE:xml_parser_benchmark>)
