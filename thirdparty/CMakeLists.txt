add_subdirectory("molfile_plugin")

# yhirose's peglib
target_include_directories(molpp SYSTEM
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp-peglib>
)
