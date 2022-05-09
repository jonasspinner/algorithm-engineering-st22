# For simplicity, this target is set up to include all C++ files placed
# in this directory. If you want to add your own executables, you'll
# need to list the appropriate sources explicitly in add_library().
# Otherwise, the benchmarks will probably complain about multiple
# definitions of main().

file(GLOB_RECURSE srcs CONFIGURE_DEPENDS
        *.cpp *.cxx *.c++ *.cc *.C *.CPP)

if (NOT srcs)
    add_library(implementation INTERFACE)
    target_include_directories(implementation INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(implementation INTERFACE Threads::Threads)
else ()
    file(GLOB_RECURSE headers CONFIGURE_DEPENDS
            *.h *.hpp *.inl *.hxx *.h++ *.hh *.H *.HPP)

    # The benchmark executables need the parameters included here,
    # so don't remove them if you list your sources manually.
    add_library(implementation STATIC
            ${srcs} ${headers}
            mst_construction_parameters.hpp)

    set_target_properties(implementation PROPERTIES LINKER_LANGUAGE CXX)
    target_include_directories(implementation INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(implementation PUBLIC Threads::Threads)
endif ()  
