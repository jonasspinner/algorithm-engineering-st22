if(CMAKE_HOST_UNIX)
  add_library(imported_binaries STATIC IMPORTED)
  set_target_properties(imported_binaries PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/binaries/libbinaries.a)
elseif(CMAKE_HOST_WIN32)
    message("Host system is windows -> no support yet")
elseif(CMAKE_HOST_APPLE)
  message("Host system is macOS -> no support yet")
endif()

add_executable(mst_construction_benchmark mst_construction_benchmark.cpp)
target_link_libraries(mst_construction_benchmark PRIVATE base_target implementation imported_binaries)

add_executable(mst_verification_benchmark mst_verification_benchmark.cpp)
target_link_libraries(mst_verification_benchmark PRIVATE base_target implementation imported_binaries)
