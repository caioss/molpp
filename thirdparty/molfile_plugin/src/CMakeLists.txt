set(MOLFILE_SOURCES
    pdbplugin.c
    mol2plugin.C
)

foreach(FILE ${MOLFILE_SOURCES})
    get_filename_component(NAME ${FILE} NAME_WLE)
    set_property(SOURCE ${FILE} PROPERTY
        COMPILE_DEFINITIONS
            VMDPLUGIN=${NAME}
    )
endforeach(FILE)

# Library defined here due to scope of set_property
add_library(molfile_plugin STATIC
    ${MOLFILE_SOURCES}
)
