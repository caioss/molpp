target_sources(molpp
    PRIVATE
        FILE_SET
            molpp_private_headers
        FILES
            DSSP.hpp
            structure.hpp

    PRIVATE
        DSSP.cpp
        structure.cpp
)
