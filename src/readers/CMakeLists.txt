target_sources(molpp
    PRIVATE
        FILE_SET
            molpp_private_headers
        FILES
            MolReader.hpp
            molfile.h
            ResidueDetect.hpp
            MolfileReader.hpp

    PRIVATE
        MolReader.cpp
        MolfileReader.cpp
)
