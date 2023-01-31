#ifndef MOLFILE_H
#define MOLFILE_H

/*
 * Molfile declarations
 */

#include "vmdplugin.h"

VMDPLUGIN_EXTERN int pdbplugin_init();
VMDPLUGIN_EXTERN int pdbplugin_register(void *v, vmdplugin_register_cb cb);
VMDPLUGIN_EXTERN int mol2plugin_init();
VMDPLUGIN_EXTERN int mol2plugin_register(void *v, vmdplugin_register_cb cb);
VMDPLUGIN_EXTERN int psfplugin_init();
VMDPLUGIN_EXTERN int psfplugin_register(void *v, vmdplugin_register_cb cb);
VMDPLUGIN_EXTERN int gromacsplugin_init();
VMDPLUGIN_EXTERN int gromacsplugin_register(void *v, vmdplugin_register_cb cb);

#endif // MOLFILE_H
