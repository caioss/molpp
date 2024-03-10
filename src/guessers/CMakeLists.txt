target_sources(molpp
    PRIVATE
        FILE_SET
            molpp_private_headers
        FILES
            ResidueBondGuesser.hpp
            AtomBondGuesser.hpp

    PRIVATE
        ResidueBondGuesser.cpp
        AtomBondGuesser.cpp
)
