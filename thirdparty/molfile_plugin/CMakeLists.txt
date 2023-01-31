target_include_directories(molpp SYSTEM
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Used by the main CMakeLists file.
set(MOLFILE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pdbplugin.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mol2plugin.C
    ${CMAKE_CURRENT_SOURCE_DIR}/src/psfplugin.c
    CACHE INTERNAL ""
)

target_sources(molpp PRIVATE
    ${MOLFILE_SOURCES}
)
