target_sources(molpp
    PRIVATE
        FILE_SET
            molpp_private_headers
        FILES
            ResiduesTable.hpp

    PRIVATE
        ElementsTable.cpp
        ResiduesTable.cpp
)
