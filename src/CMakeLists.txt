target_sources(molpp
    PRIVATE
        FILE_SET
            molpp_private_headers
        TYPE
            HEADERS
)

add_subdirectory("core")
add_subdirectory("selections")
add_subdirectory("readers")
add_subdirectory("tables")
add_subdirectory("guessers")
add_subdirectory("tools")
add_subdirectory("analysis")
