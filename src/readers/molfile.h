#ifndef MOLFILE_H
#define MOLFILE_H

/*
 * Molfile declarations
 */

#include "vmdplugin.h"

VMDPLUGIN_EXTERN int pdbplugin_init();
VMDPLUGIN_EXTERN int pdbplugin_register(void *v, vmdplugin_register_cb cb);

#endif // MOLFILE_H
