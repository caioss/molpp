/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2016 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: jsplugin.c,v $
 *      $Author: johns $       $Locker:  $             $State: Exp $
 *      $Revision: 1.91 $       $Date: 2020/10/21 16:01:57 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   This plugin implements a high-performance binary molecular structure
 *   and trajectory storage format.  This file format currently uses a simple 
 *   non-redundant hash table approach for compression of per-atom
 *   character strings, properties, and tags.  Trajectory data is stored 
 *   with a file structure that avoids the need to transpose or convert
 *   dense blocks of cartesian coordinates into the most commonly used
 *   interleaved  x/y/z ordering.  The file structure also enables zero-copy 
 *   vectorized I/O methods to be used high-performance reads for 
 *   visualization and analysis with reduced operating system overhead.
 *   The plugin optionally supports the use of a block-based file structure
 *   and block-aligned memory buffers for direct I/O that bypasses the 
 *   OS filesystem buffer caches for multi-gigabyte-per-second read rates
 *   from SSD RAID arrays.
 *
 *   At present, only VMD, NAMD, and psfgen make use of this format.
 *   It started out as a test/example code and is slowly becoming
 *   more robust and featureful.
 *
 *   We should be able to implement a selective read approach that gathers
 *   discontiguous blocks of the trajectory using the POSIX lio_listio()
 *   APIs.  On Unix we currently use I/O wrappers that are based on the 
 *   lseek() and readv() APIs.  By using lio_listio() we could eliminate
 *   the separate lseek calls and service multiple timestep reads in a 
 *   single request, even included cases with discontiguous requests.
 *
 *   VMD test results for Linux host with an 8-way RAID0 of commodity 
 *   Intel 510 SSDs with SATA III 6Gbit/sec interfaces:
 *     Non-direct I/O using standard readv(): 1203 MB/sec
 *     Direct I/O, readv(), 4KB blocked file: 2130 MB/sec
 ***************************************************************************/
/***
 ** Standalone test binary compilation flags for 64-bit Linux:
    cc -O3 -m64 -I../../include -DTEST_JSPLUGIN jsplugin.c \
      -o ~/bin/readjs -lm
 **
 ** Standalone test binary compilation flags for 64-bit Linux w/ CUDA:
    cc -O3 -m64 -I../../include -I/usr/local/cuda/include \
      -DTEST_JSPLUGIN -DENABLECUDATESTS jsplugin.c \
      -o ~/bin/readjs -L/usr/local/cuda/lib64 -lcudart -lm
 **
 ** Standalone test binary compilation flags for 64-bit Linux w/ CUDA+GDS:
    cc -O3 -m64 -I../../include -I/usr/local/cuda/include \
      -I/Projects/vmd/cuda/johns/gds/gds-alpha/lib \
      -DTEST_JSPLUGIN -DENABLECUDATESTS -DENABLECUDAGDS jsplugin.c \
      -o ~/bin/readjs -L/Projects/vmd/cuda/johns/gds/gds-alpha/lib -lcufile \
      -L/usr/local/cuda/lib64 -lcudart -lm
 **
 ** Standalone test binary compilation flags for Solaris:
    cc -fast -xarch=v9a -I../../include -DTEST_JSPLUGIN jsplugin.c \
      -o ~/bin/readjs -lm
 **
 ** Profiling flags for Solaris:
    cc -xpg -fast -xarch=v9a -g -I../../include -DTEST_JSPLUGIN jsplugin.c \
      -o ~/bin/readjs -lm
 **
 ** Test run for DGX-2:
    ~/bin/readjs /raid/scratch/solvatebar1204kb.js
 ***************************************************************************/

#define INFOMSGS  1

#if 1
#define ENABLEJSSHORTREADS 1
#endif

#define VMDPLUGIN_STATIC
#include "largefiles.h"   /* platform dependent 64-bit file I/O defines */
#include "fastio.h"       /* must come before others, for O_DIRECT...   */

#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "hash.h"
#include "endianswap.h"
#include "molfile_plugin.h"


/* allocate memory and return a pointer that is aligned on a given   */
/* byte boundary, to be used for page- or sector-aligned I/O buffers */
/* We use this if posix_memalign() is not available...               */
#if defined(_WIN64)  /* sizeof(size_t) == sizeof(void*) */
#define myintptrtype size_t
#elif 1 /* sizeof(unsigned long) == sizeof(void*) */
#define myintptrtype unsigned long
#else
#define myintptrtype uintptr_t  /* C99 */
#endif
/*
 * XXX On MSVC we get warnings about type conversions for 
 *     size_t vs. fio_size_t
 */
static void *alloc_aligned_ptr(size_t sz, size_t blocksz, void **unalignedptr) {
  /* pad the allocation to an even multiple of the block size */
  size_t padsz = (sz + (blocksz - 1)) & (~(blocksz - 1));
  void * ptr = malloc(padsz + blocksz);
  *unalignedptr = ptr;
  return (void *) ((((myintptrtype) ptr) + (blocksz-1)) & (~(blocksz-1)));
}


#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define JSHEADERSTRING   "JS Binary Structure and Trajectory File Format"                
#define JSMAGICNUMBER    0x31337
#define JSENDIANISM      0x12345678

#define JSMAJORVERSION   2
#define JSMINORVERSION   19

#define JSNFRAMESOFFSET  (strlen(JSHEADERSTRING) + 20)

#define JSNOERR             0
#define JSBADFILE           1
#define JSBADFORMAT         2


/* Threshold atom count beyond which block-based I/O is used by default */
/* The overhead from block-alignment padding bytes becomes essentially  */
/* inconsequential (< 1%) for structures with more than 50,000 atoms.   */
#define JSBLOCKIO_THRESH    50000


/* Option flag macros and their meanings */
#define JSOPT_NOOPTIONS     0x00000000  /* no structure, only coords    */

/* Timesteps are block-size padded and page- or sector-aligned for      */
/* direct I/O, using  OS-specific APIs that completely bypass the OS    */
/* kernel filesystem buffer cache.                                      */
/* The use of direct I/O APIs can raise performance tremendously on     */
/* high-end RAIDs.  Tests on an 8-way RAID0 of Intel 510 SSDs raise the */
/* peak I/O rate from 1100 MB/sec up to 2020 MB/sec with direct I/O.    */
#define JSOPT_TS_BLOCKIO    0x10000000

/* large data blocks */
#define JSOPT_STRUCTURE     0x00000001  /* file contains structure data */
#define JSOPT_BONDS         0x00000002  /* file contains bond info      */
#define JSOPT_BONDORDERS    0x00000004  /* file contains bond orders    */
#define JSOPT_ANGLES        0x00000008  /* file contains angle info     */
#define JSOPT_CTERMS        0x00000010  /* file contains cross-terms    */

/* optional per-atom fields */
#define JSOPT_OCCUPANCY     0x00000100  /* file contains occupancy      */
#define JSOPT_BFACTOR       0x00000200  /* file contains b-factor       */
#define JSOPT_MASS          0x00000400  /* file contains masses         */
#define JSOPT_CHARGE        0x00000800  /* file contains charges        */
#define JSOPT_RADIUS        0x00001000  /* file contains radii          */
#define JSOPT_ATOMICNUMBER  0x00002000  /* file contains atomic numbers */

typedef struct {
  int verbose;                 /* flag to enable console info output    */
  fio_fd fd;                   /* main file descriptor                  */
  ptrdiff_t natoms;            /* handle uses a long type for natoms to */
                               /* help force promotion of file offset   */
                               /* arithmetic to long types              */

#if JSMAJORVERSION > 1
  int parsed_structure;        /* flag indicating structure is parsed   */
  char *path;                  /* path to file                          */

  /* info for block-based direct I/O */ 
  int directio_pgsize_queried; /* caller has queried page/blocksz       */
  int directio_enabled;        /* block-based direct I/O is available   */
  fio_fd directio_fd;          /* block-based direct I/O using O_DIRECT */
  int directio_block_size;     /* block size to use for direct ts I/O   */
  void *directio_ucell_ptr;    /* unaligned unit cell buffer ptr        */
  void *directio_ucell_blkbuf; /* block-aligned unit cell buffer pt r   */

  /* timestep file offset, block padding, and stride information */
  fio_size_t ts_file_offset;   /* file offset to first timestep         */
  fio_size_t ts_crd_sz;        /* size of TS coordinates                */
  fio_size_t ts_crd_padsz;     /* size of TS block-padded coordinates   */
  fio_size_t ts_ucell_sz;      /* size of TS unit cell                  */
  fio_size_t ts_ucell_padsz;   /* size of TS block-padded unit cell     */
  
  /* structure info */
  int optflags;
  molfile_atom_t *atomlist;
  molfile_metadata_t *meta;

  /* bond info */
  int nbonds;
  int *bondfrom;
  int *bondto;
  float *bondorders;

  /* angle/dihedral/improper/cross-term info */
  int numangles, *angles;
  int numdihedrals, *dihedrals;
  int numimpropers, *impropers;
  int numcterms, *cterms;
#endif

  /* trajectory info */
  int nframes;
  double tsdelta;
  int reverseendian;
  int with_unitcell;

  /* convenient buffer full of zeros for block-multiple padding */
  unsigned char blockpad[MOLFILE_DIRECTIO_MAX_BLOCK_SIZE];
} jshandle;


/* report the block size required to read this JS file */
static int read_js_timestep_pagealign_size(void *v, int *pagealignsz) {
  jshandle *js = (jshandle *)v;

  // mark that the caller has queried the page alignment size
  js->directio_pgsize_queried = 1;

  /* assigne page alignment size based on file contents */
  if (js->optflags & JSOPT_TS_BLOCKIO) 
    *pagealignsz = js->directio_block_size;
  else 
    *pagealignsz = 1;

  return 0;
}


/* Use block-based I/O by default when writing structures larger */
/* than JSBLOCKIO_THRESH atoms, or when directed by the user and */
/* not otherwise prohibited...                                   */
static void js_blockio_check_and_set(jshandle *js) {
  if ((getenv("VMDJSNOBLOCKIO") == NULL) && 
      ((js->natoms > JSBLOCKIO_THRESH) || getenv("VMDJSBLOCKIO"))) {
    js->optflags |= JSOPT_TS_BLOCKIO;
    js->directio_block_size = MOLFILE_DIRECTIO_MIN_BLOCK_SIZE; 
  }
}


static void *open_js_read(const char *path, const char *filetype, int *natoms) {
  jshandle *js;
  int jsmagicnumber, jsendianism, jsmajorversion, jsminorversion;
  struct stat stbuf;
  char strbuf[1024];
  int tmpnatoms=0;

  if (!path) return NULL;

  /* See if the file exists, and get its size */
  memset(&stbuf, 0, sizeof(struct stat));
  if (stat(path, &stbuf)) {
    printf("jsplugin) Could not access file '%s'.\n", path);
    perror("jsplugin) stat: ");
/*    return NULL; */
  }

  js = (jshandle *)malloc(sizeof(jshandle));
  memset(js, 0, sizeof(jshandle));
  js->verbose = (getenv("VMDJSVERBOSE") != NULL);
#if defined(_WIN64)
  js->verbose = 1;
#endif

#if JSMAJORVERSION > 1
  js->parsed_structure=0;
  js->directio_block_size=1;
  js->directio_ucell_ptr = NULL;
  js->directio_ucell_blkbuf = NULL;

  js->directio_pgsize_queried=0;
  js->directio_enabled=0;
  js->ts_file_offset=0;
  js->ts_crd_sz=0;
  js->ts_ucell_sz=0;
  js->ts_crd_padsz=0;
  js->ts_ucell_padsz=0;
#endif

  if (fio_open(path, FIO_READ, &js->fd) < 0) {
    printf("jsplugin) Could not open file '%s' for reading.\n", path);
    free(js);
    return NULL;
  }

  /* emit header information */
  fio_fread(strbuf, strlen(JSHEADERSTRING), 1, js->fd);
  strbuf[strlen(JSHEADERSTRING)] = '\0';
  if (strcmp(strbuf, JSHEADERSTRING)) {
    printf("jsplugin) Bad trajectory header!\n");
    printf("jsplugin) Read string: %s\n", strbuf);
    fio_fclose(js->fd);
    free(js);
    return NULL;
  }

  fio_read_int32(js->fd, &jsmagicnumber);
  fio_read_int32(js->fd, &jsendianism);
  fio_read_int32(js->fd, &jsmajorversion);
  fio_read_int32(js->fd, &jsminorversion);
  fio_read_int32(js->fd, &tmpnatoms); /* handle-internal natoms is a long */
  fio_read_int32(js->fd, &js->nframes);
  if ((jsmagicnumber != JSMAGICNUMBER) || (jsendianism != JSENDIANISM)) {
#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin) opposite endianism file, enabling byte swapping\n");
#endif
    js->reverseendian = 1;
    swap4_aligned(&jsmagicnumber, 1);
    swap4_aligned(&jsendianism, 1);
    swap4_aligned(&jsmajorversion, 1);
    swap4_aligned(&jsminorversion, 1);
    swap4_aligned(&tmpnatoms, 1);
    swap4_aligned(&js->nframes, 1);
  } else {
#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin) native endianism file\n");
#endif
  }

  if ((jsmagicnumber != JSMAGICNUMBER) || (jsendianism != JSENDIANISM)) {
    fio_fclose(js->fd);
    free(js);
    return NULL;
  }
 
  if (jsmajorversion != JSMAJORVERSION) {
    printf("jsplugin) major version mismatch\n");
    printf("jsplugin)   file version: %d\n", jsmajorversion);
    printf("jsplugin)   plugin version: %d\n", JSMAJORVERSION);
    fio_fclose(js->fd);
    free(js);
    return NULL;
  }

  /* Copy integer natoms to handle natoms, could be a long. */
  /* The handle natoms uses long to help force promotion of */
  /* integer file offset calculations to long types...      */
  js->natoms = tmpnatoms;
  *natoms = tmpnatoms;

  /* copy path if we succeeded in opening the file */
  js->path = (char *) calloc(strlen(path)+1, 1);
  strcpy(js->path, path);

#if 1
  /* read flags data from the file */
  fio_read_int32(js->fd, &js->optflags); 
  if (js->reverseendian)
    swap4_aligned(&js->optflags, 1);

#if defined(INFOMSGS)
  if (js->verbose)
    printf("jsplugin) read option flags: %0x08x\n", js->optflags);
#endif

  /* Check to see if block-based trajectory I/O is used  */
  /* and read in the block size for this file.           */
  if (js->optflags & JSOPT_TS_BLOCKIO) {
    fio_fread(&js->directio_block_size, sizeof(int), 1, js->fd);
    if (js->reverseendian)
      swap4_aligned(&js->directio_block_size, 1);

#if defined(INFOMSGS)
    if (js->verbose) {
      printf("jsplugin) File uses direct I/O block size: %d bytes\n", 
             js->directio_block_size);
    }
#endif

    /* Check to ensure that we can handle the block size used by the */
    /* file we are reading.  We may use variable block sizes in      */
    /* the future as more high-end filesystems begin to support      */
    /* 8KB, 16KB, or larger block sizes for enhanced sequential I/O  */
    /* performance on very large files.                              */
    if (js->directio_block_size > MOLFILE_DIRECTIO_MAX_BLOCK_SIZE) {
      printf("jsplugin) File block size exceeds jsplugin block size limit.\n");
      printf("jsplugin) Direct I/O unavailable for file '%s'\n", js->path);
    } else {
      if (fio_open(js->path, FIO_READ | FIO_DIRECT, &js->directio_fd) < 0) {
        printf("jsplugin) Direct I/O unavailable for file '%s'\n", js->path);
      } else {
        js->directio_enabled = 1;
      } 
    }
  }

#if defined(ENABLEJSSHORTREADS)
  /* test code for an implementation that does short reads that */
  /* skip bulk solvent, useful for faster loading of very large */
  /* structures                                                 */
  if (getenv("VMDJSMAXATOMIDX") != NULL) {
    ptrdiff_t maxatomidx = atoi(getenv("VMDJSMAXATOMIDX"));
    if (maxatomidx < 0)
      maxatomidx = 0;
    if (maxatomidx >= js->natoms)
      maxatomidx = js->natoms - 1;

    printf("jsplugin) Short-reads of timesteps enabled: %td / %td atoms (%.2f%%)\n",
           maxatomidx, js->natoms, 100.0*(maxatomidx+1) / ((double) js->natoms));
  }
#endif
#endif

  return js;
}


#if JSMAJORVERSION > 1

/* Compute the file offset for the first timestep and move */
/* the file pointer to the correct position to read/write  */
/* the first timestep.  Takes care of block alignment when */
/* needed.                                                 */ 
static int js_calc_timestep_blocking_info(void *mydata) {
  fio_size_t ts_block_offset, bszmask;
  jshandle *js = (jshandle *) mydata;
  int iorc=0;

  /* Record the current file offset so we can use it to */
  /* compute the absolute offset to the first timestep. */
  js->ts_file_offset = fio_ftell(js->fd);

  /* pad current offset to the start of the next block  */ 
  bszmask = js->directio_block_size - 1;
  ts_block_offset = (js->ts_file_offset + bszmask) & (~bszmask);

#if defined(INFOMSGS)
  if (js->verbose) {
    printf("jsplugin) TS block size %td  curpos: %td  blockpos: %td\n", 
           (ptrdiff_t) js->directio_block_size, 
           (ptrdiff_t) js->ts_file_offset, 
           (ptrdiff_t) ts_block_offset);
  }
#endif

  /* seek to the first block of the first timestep */
  js->ts_file_offset = ts_block_offset;
  if (js->directio_enabled)
    iorc = fio_fseek(js->directio_fd, js->ts_file_offset, FIO_SEEK_SET);
  else
    iorc = fio_fseek(js->fd, js->ts_file_offset, FIO_SEEK_SET);
  if (iorc < 0) {
    perror("jsplugin) fseek(): ");
  }

  /* compute timestep block padding/skipping for both */
  /* coordinate blocks and unit cell blocks           */
  js->ts_crd_sz = js->natoms * 3L * sizeof(float);
  js->ts_crd_padsz = (js->ts_crd_sz + bszmask) & (~bszmask);

  js->ts_ucell_sz = 6L * sizeof(double);
  js->ts_ucell_padsz = (js->ts_ucell_sz + bszmask) & (~bszmask);

  /* allocate TS unit cell buffer in an aligned, block-size-multiple buffer */
  /* unaligned unit cell buffer ptr */
#if defined(USE_POSIX_MEMALIGN)
  if (posix_memalign((void**) &js->directio_ucell_ptr, 
      js->directio_block_size, js->ts_ucell_padsz)) {
    printf("jsplugin) Couldn't allocate aligned unit cell block buffer!\n");
  }
  /* the returned pointer is already block-aligned, and can free() */
  js->directio_ucell_blkbuf = js->directio_ucell_ptr;
#else
  js->directio_ucell_blkbuf = (float *) 
    alloc_aligned_ptr(js->ts_ucell_padsz, js->directio_block_size, 
                      (void**) &js->directio_ucell_ptr);
#endif

#if defined(INFOMSGS)
  if (js->verbose) {
    printf("jsplugin) TS crds sz: %td psz: %td  ucell sz: %td psz: %td\n",
           (ptrdiff_t) js->ts_crd_sz,
           (ptrdiff_t) js->ts_crd_padsz, 
           (ptrdiff_t) js->ts_ucell_sz, 
           (ptrdiff_t) js->ts_ucell_padsz);
  }
#endif

  return MOLFILE_SUCCESS;
}


static int read_js_structure(void *mydata, int *optflags,
                             molfile_atom_t *atoms) {
  jshandle *js = (jshandle *) mydata;
  ptrdiff_t i;

  if (optflags != NULL)
    *optflags = MOLFILE_NOOPTIONS; /* set to no options until we read them */

#if 0
  /* read flags data from the file */
  fio_read_int32(js->fd, &js->optflags); 
  if (js->reverseendian)
    swap4_aligned(&js->optflags, 1);

#if defined(INFOMSGS)
  if (js->verbose)
    printf("jsplugin) read option flags: %0x08x\n", js->optflags);
#endif

  /* Check to see if block-based trajectory I/O is used  */
  /* and read in the block size for this file.           */
  if (js->optflags & JSOPT_TS_BLOCKIO) {
    fio_fread(&js->directio_block_size, sizeof(int), 1, js->fd);
    if (js->reverseendian)
      swap4_aligned(&js->directio_block_size, 1);

#if defined(INFOMSGS)
    if (js->verbose) {
      printf("jsplugin) File uses direct I/O block size: %d bytes\n", 
             js->directio_block_size);
    }
#endif

    /* Check to ensure that we can handle the block size used by the */
    /* file we are reading.  We may use variable block sizes in      */
    /* the future as more high-end filesystems begin to support      */
    /* 8KB, 16KB, or larger block sizes for enhanced sequential I/O  */
    /* performance on very large files.                              */
    if (js->directio_block_size > MOLFILE_DIRECTIO_MAX_BLOCK_SIZE) {
      printf("jsplugin) File block size exceeds jsplugin block size limit.\n");
      printf("jsplugin) Direct I/O unavailable for file '%s'\n", js->path);
    } else {
      if (fio_open(js->path, FIO_READ | FIO_DIRECT, &js->directio_fd) < 0) {
        printf("jsplugin) Direct I/O unavailable for file '%s'\n", js->path);
      } else {
        js->directio_enabled = 1;
      } 
    }
  }
#endif


  /* emit warning message if the caller didn't check the required */
  /* alignment size, but the file supports block based direct I/O */
  if (js->directio_enabled && !js->directio_pgsize_queried) {
    printf("jsplugin) Warning: File supports block-based direct I/O, but\n");
    printf("jsplugin)          caller failed to query required alignment.\n");
    printf("jsplugin)          Block-based direct I/O is now disabled.\n");
     
    js->directio_enabled=0; // ensure we disable direct I/O early on
  }

#if defined(INFOMSGS)
  if (js->verbose) {
    printf("jsplugin) Direct I/O %sabled for file '%s'\n", 
           (js->directio_enabled) ? "en" : "dis", js->path);
  }
#endif


#if 0
#if defined(ENABLEJSSHORTREADS)
  /* test code for an implementation that does short reads that */
  /* skip bulk solvent, useful for faster loading of very large */
  /* structures                                                 */
  if (getenv("VMDJSMAXATOMIDX") != NULL) {
    ptrdiff_t maxatomidx = atoi(getenv("VMDJSMAXATOMIDX"));
    if (maxatomidx < 0)
      maxatomidx = 0;
    if (maxatomidx >= js->natoms)
      maxatomidx = js->natoms - 1;

    printf("jsplugin) Short-reads of timesteps enabled: %ld / %ld atoms (%.2f%%)\n",
           maxatomidx, js->natoms, 100.0*(maxatomidx+1) / ((double) js->natoms));
  }
#endif
#endif

  /* Mark the handle to indicate we've parsed the structure.             */
  /* If any errors occur after this point, they are likely fatal anyway. */
  js->parsed_structure = 1;

  /* determine whether or not this file contains structure info or not */
  if (js->optflags & JSOPT_STRUCTURE) {
    int numatomnames, numatomtypes, numresnames, numsegids, numchains;
    char **atomnames = NULL;
    char **atomtypes = NULL;
    char **resnames = NULL;
    char **segids = NULL;
    char **chains = NULL;
    short *shortbuf = NULL; /* temp buf for decoding atom records */
    int *intbuf = NULL;     /* temp buf for decoding atom records */
    float *fltbuf = NULL;   /* temp buf for decoding atom records */
 
    /* read in block of name string table sizes */
    fio_read_int32(js->fd, &numatomnames); 
    fio_read_int32(js->fd, &numatomtypes); 
    fio_read_int32(js->fd, &numresnames);
    fio_read_int32(js->fd, &numsegids);
    fio_read_int32(js->fd, &numchains); 
    if (js->reverseendian) {
      swap4_aligned(&numatomnames, 1);
      swap4_aligned(&numatomtypes, 1);
      swap4_aligned(&numresnames, 1);
      swap4_aligned(&numsegids, 1);
      swap4_aligned(&numchains, 1);
    }

#if defined(INFOMSGS)
    if (js->verbose) {
      printf("jsplugin) reading string tables...\n");
      printf("jsplugin) %d %d %d %d %d\n",
             numatomnames, numatomtypes, numresnames, numsegids, numchains);
    }
#endif

    /* skip forward to first TS if the caller gives us NULL ptrs */
    if (optflags == NULL && atoms == NULL) {
      size_t offset=0;
      offset += numatomnames * (16L * sizeof(char));
      offset += numatomtypes * (16L * sizeof(char));
      offset += numresnames  * (8L * sizeof(char));
      offset += numsegids    * (8L * sizeof(char));
      offset += numchains    * (2L * sizeof(char));
      offset += js->natoms * sizeof(short); /* atom name indices    */
      offset += js->natoms * sizeof(short); /* atom type indices    */
      offset += js->natoms * sizeof(short); /* residue name indices */
      offset += js->natoms * sizeof(short); /* segment name indices */
      offset += js->natoms * sizeof(short); /* chain name indices   */
      offset += js->natoms * sizeof(int);   /* residue indices      */
      
      /* optional per-atom fields */
      if (js->optflags & JSOPT_OCCUPANCY)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_BFACTOR)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_MASS)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_CHARGE)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_RADIUS)
        offset += js->natoms * sizeof(float); 
      if (js->optflags & JSOPT_ATOMICNUMBER)
        offset += js->natoms * sizeof(int);

      fio_fseek(js->fd, offset, FIO_SEEK_CUR);
      offset=0;

      /* these require actually seeking as we process... */
      if (js->optflags & JSOPT_BONDS) {
        fio_fread(&js->nbonds, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->nbonds, 1);
#if defined(INFOMSGS)
        if (js->verbose) {
          printf("jsplugin)   %d bonds...\n", js->nbonds);
        }
#endif

        offset += 2L * js->nbonds * sizeof(int);
        if (js->optflags & JSOPT_BONDORDERS)
          offset += js->nbonds * sizeof(float);

        fio_fseek(js->fd, offset, FIO_SEEK_CUR);
        offset=0;
      }

      if (js->optflags & JSOPT_ANGLES) {
        fio_fread(&js->numangles, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->numangles, 1);
#if defined(INFOMSGS)
        if (js->verbose) {
          printf("jsplugin)   %d angles...\n", js->numangles);
        }
#endif
        fio_fseek(js->fd, sizeof(int)*3L*js->numangles, FIO_SEEK_CUR);

        fio_fread(&js->numdihedrals, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->numdihedrals, 1);
#if defined(INFOMSGS)
        if (js->verbose) {
          printf("jsplugin)   %d dihedrals...\n", js->numdihedrals);
        }
#endif
        fio_fseek(js->fd, sizeof(int)*4L*js->numdihedrals, FIO_SEEK_CUR);

        fio_fread(&js->numimpropers, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->numimpropers, 1);
#if defined(INFOMSGS)
        if (js->verbose) {
          printf("jsplugin)   %d impropers...\n", js->numimpropers);
        }
#endif
        fio_fseek(js->fd, sizeof(int)*4L*js->numimpropers, FIO_SEEK_CUR);
      }

      if (js->optflags & JSOPT_CTERMS) {
        fio_fread(&js->numcterms, sizeof(int), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(&js->numcterms, 1);
#if defined(INFOMSGS)
        if (js->verbose) {
          printf("jsplugin)   %d cterms...\n", js->numcterms);
        }
#endif
        fio_fseek(js->fd, sizeof(int)*8L*js->numcterms, FIO_SEEK_CUR);
      }
  
      /* record the file offset for the first timestep */
      js_calc_timestep_blocking_info(js);

      return MOLFILE_SUCCESS;
    }


    /* allocate string tables */
    atomnames = (char **) malloc(numatomnames * sizeof(char *));
    atomtypes = (char **) malloc(numatomtypes * sizeof(char *));
    resnames  = (char **) malloc(numresnames  * sizeof(char *));
    segids    = (char **) malloc(numsegids    * sizeof(char *));
    chains    = (char **) malloc(numchains    * sizeof(char *));

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   atom names...\n");
#endif

    /* read in the string tables */
    for (i=0; i<numatomnames; i++) {
      atomnames[i] = (char *) malloc(16L * sizeof(char));
      fio_fread(atomnames[i], 16L * sizeof(char), 1, js->fd);
    }

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   atom types...\n");
#endif
    for (i=0; i<numatomtypes; i++) {
      atomtypes[i] = (char *) malloc(16L * sizeof(char));
      fio_fread(atomtypes[i], 16L * sizeof(char), 1, js->fd);
    }

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   residue names...\n");
#endif
    for (i=0; i<numresnames; i++) {
      resnames[i] = (char *) malloc(8L * sizeof(char));
      fio_fread(resnames[i], 8L * sizeof(char), 1, js->fd);
    }

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   segment names...\n");
#endif
    for (i=0; i<numsegids; i++) {
      segids[i] = (char *) malloc(8L * sizeof(char));
      fio_fread(segids[i], 8L * sizeof(char), 1, js->fd);
    }

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   chain names...\n");
#endif
    for (i=0; i<numchains; i++) {
      chains[i] = (char *) malloc(2L * sizeof(char));
      fio_fread(chains[i], 2L * sizeof(char), 1, js->fd);
    }

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin) reading numeric field tables...\n");
#endif
    /* read in all of the atom fields */
    shortbuf = (short *) malloc(js->natoms * sizeof(short));

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   atom name indices...\n");
#endif
    /* read in atom names */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].name, atomnames[shortbuf[i]]);
    }
    for (i=0; i<numatomnames; i++)
      free(atomnames[i]);
    free(atomnames);

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   atom type indices...\n");
#endif
    /* read in atom types */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].type, atomtypes[shortbuf[i]]);
    }
    for (i=0; i<numatomtypes; i++)
      free(atomtypes[i]);
    free(atomtypes);

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   residue name indices...\n");
#endif
    /* read in resnames */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].resname, resnames[shortbuf[i]]);
    }
    for (i=0; i<numresnames; i++)
      free(resnames[i]);
    free(resnames);
    
#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   segment name indices...\n");
#endif
    /* read in segids */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].segid, segids[shortbuf[i]]);
    }
    for (i=0; i<numsegids; i++)
      free(segids[i]);
    free(segids);

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   chain name indices...\n");
#endif
    /* read in chains */
    fio_fread(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    if (js->reverseendian)
      swap2_aligned(shortbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      strcpy(atoms[i].chain, chains[shortbuf[i]]);
    }
    for (i=0; i<numchains; i++)
      free(chains[i]);
    free(chains);

    if (shortbuf != NULL) {
      free(shortbuf);
      shortbuf=NULL;
    }

    /* 
     * read in integer data blocks 
     */
    intbuf = (int *) malloc(js->natoms * sizeof(int));

#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin)   residue indices...\n");
#endif
    /* read in resid */
    fio_fread(intbuf, js->natoms * sizeof(int), 1, js->fd);
    if (js->reverseendian)
      swap4_aligned(intbuf, js->natoms);
    for (i=0; i<js->natoms; i++) {
      atoms[i].resid = intbuf[i];
    }    
     
    if (intbuf != NULL) {
      free(intbuf);
      intbuf = NULL;
    }


#if defined(INFOMSGS)
    if (js->verbose)
      printf("jsplugin) reading optional per-atom tables...\n");
#endif
    /*
     * read in optional single-precision float data blocks
     */ 
    if (js->optflags & (JSOPT_OCCUPANCY | JSOPT_BFACTOR | 
        JSOPT_MASS | JSOPT_RADIUS | JSOPT_CHARGE)) 
      fltbuf = (float *) malloc(js->natoms * sizeof(float));

    /* read in optional data if it exists */
    if (js->optflags & JSOPT_OCCUPANCY) {
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   occupancy...\n");
#endif
      *optflags |= MOLFILE_OCCUPANCY;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].occupancy = fltbuf[i];
      }    
    }

    if (js->optflags & JSOPT_BFACTOR) {
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   bfactor...\n");
#endif
      *optflags |= MOLFILE_BFACTOR;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].bfactor = fltbuf[i];
      }    
    }

    if (js->optflags & JSOPT_MASS) { 
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   mass...\n");
#endif
      *optflags |= MOLFILE_MASS;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].mass = fltbuf[i];
      }    
    }

    if (js->optflags & JSOPT_CHARGE) { 
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   charge...\n");
#endif
      *optflags |= MOLFILE_CHARGE;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].charge = fltbuf[i];
      }    
    }

    if (js->optflags & JSOPT_RADIUS) { 
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   radius...\n");
#endif
      *optflags |= MOLFILE_RADIUS;
      fio_fread(fltbuf, js->natoms * sizeof(float), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(fltbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].radius = fltbuf[i];
      }    
    }

    if (fltbuf != NULL) {
      free(fltbuf);
      fltbuf=NULL;
    }

    /*
     * read in optional integer data blocks
     */ 
    if (js->optflags & JSOPT_ATOMICNUMBER)
      intbuf = (int *) malloc(js->natoms * sizeof(int));

    if (js->optflags & JSOPT_ATOMICNUMBER) { 
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   atomic number...\n");
#endif
      *optflags |= MOLFILE_ATOMICNUMBER;
      fio_fread(intbuf, js->natoms * sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(intbuf, js->natoms);
      for (i=0; i<js->natoms; i++) {
        atoms[i].atomicnumber = intbuf[i];
      }    
    }

    if (intbuf != NULL) {
      free(intbuf);
      intbuf = NULL;
    }


    /*
     * read in bonds and fractional bond orders
     */ 
    if (js->optflags & JSOPT_BONDS) {
      fio_fread(&js->nbonds, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->nbonds, 1);
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   %d bonds...\n", js->nbonds);
#endif

      js->bondfrom = (int *) malloc(js->nbonds * sizeof(int));
      js->bondto = (int *) malloc(js->nbonds * sizeof(int));
      fio_fread(js->bondfrom, js->nbonds * sizeof(int), 1, js->fd);
      fio_fread(js->bondto, js->nbonds * sizeof(int), 1, js->fd);
      if (js->reverseendian) {
        swap4_aligned(js->bondfrom, js->nbonds);
        swap4_aligned(js->bondto, js->nbonds);
      }

      if (js->optflags & JSOPT_BONDORDERS) {
#if defined(INFOMSGS)
        if (js->verbose)
          printf("jsplugin)   bond orders...\n");
#endif
        js->bondorders = (float *) malloc(js->nbonds * sizeof(float));
        fio_fread(js->bondorders, js->nbonds * sizeof(float), 1, js->fd);
        if (js->reverseendian)
          swap4_aligned(js->bondorders, js->nbonds);
      }
    }

    if (js->optflags & JSOPT_ANGLES) {
      fio_fread(&js->numangles, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numangles, 1);
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   %d angles...\n", js->numangles);
#endif
      js->angles = (int *) malloc(3L * js->numangles * sizeof(int));
      fio_fread(js->angles, sizeof(int)*3L*js->numangles, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(js->angles, 3L*js->numangles);

      fio_fread(&js->numdihedrals, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numdihedrals, 1);
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   %d dihedrals...\n", js->numdihedrals);
#endif
      js->dihedrals = (int *) malloc(4L * js->numdihedrals * sizeof(int));
      fio_fread(js->dihedrals, sizeof(int)*4L*js->numdihedrals, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(js->dihedrals, 4L*js->numdihedrals);

      fio_fread(&js->numimpropers, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numimpropers, 1);
      js->impropers = (int *) malloc(4L * js->numimpropers * sizeof(int));
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   %d impropers...\n", js->numimpropers);
#endif
      fio_fread(js->impropers, sizeof(int)*4L*js->numimpropers, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(js->impropers, 4L*js->numimpropers);
    }

    if (js->optflags & JSOPT_CTERMS) {
      fio_fread(&js->numcterms, sizeof(int), 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(&js->numcterms, 1);
      js->cterms = (int *) malloc(8L * js->numcterms * sizeof(int));
#if defined(INFOMSGS)
      if (js->verbose)
        printf("jsplugin)   %d cterms...\n", js->numcterms);
#endif
      fio_fread(js->cterms, sizeof(int)*8L*js->numcterms, 1, js->fd);
      if (js->reverseendian)
        swap4_aligned(js->cterms, 8L*js->numcterms);
    }

#if defined(INFOMSGS)
    if (js->verbose) {
      printf("jsplugin) final optflags: %08x\n", *optflags);
      printf("jsplugin) structure information complete\n");
    }
#endif

    /* record the file offset for the first timestep */
    js_calc_timestep_blocking_info(js);

    return MOLFILE_SUCCESS;
  }

#if defined(INFOMSGS)
  if (js->verbose)
    printf("jsplugin) no structure information available\n");
#endif

  /* record the file offset for the first timestep */
  js_calc_timestep_blocking_info(js);

  /* else, we have no structure information */
  return MOLFILE_NOSTRUCTUREDATA;
}


static int read_js_bonds(void *v, int *nbonds, int **fromptr, int **toptr, 
                         float **bondorder, int **bondtype, 
                         int *nbondtypes, char ***bondtypename) {
  jshandle *js = (jshandle *)v;

  *nbonds = 0;
  *fromptr = NULL;
  *toptr = NULL;
  *bondorder = NULL;
  *bondtype = NULL;
  *nbondtypes = 0;
  *bondtypename = NULL;

  if (js->optflags & JSOPT_BONDS) {
    *nbonds = js->nbonds;
    *fromptr = js->bondfrom;
    *toptr = js->bondto;

    if (js->optflags & JSOPT_BONDORDERS) {
      *bondorder = js->bondorders;
    }
  }

  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 14
static int read_js_angles(void *v, int *numangles, int **angles, 
                          int **angletypes, int *numangletypes, 
                          char ***angletypenames, int *numdihedrals,
                          int **dihedrals, int **dihedraltypes, 
                          int *numdihedraltypes, char ***dihedraltypenames,
                          int *numimpropers, int **impropers, 
                          int **impropertypes, int *numimpropertypes, 
                          char ***impropertypenames, int *numcterms, 
                          int **cterms, int *ctermcols, int *ctermrows) {
  jshandle *js = (jshandle *)v;

  /* initialize data to zero */
  *numangles         = 0;
  *angles            = NULL;
  *angletypes        = NULL;
  *numangletypes     = 0;
  *angletypenames    = NULL;
  *numdihedrals      = 0;
  *dihedrals         = NULL;
  *dihedraltypes     = NULL;
  *numdihedraltypes  = 0;
  *dihedraltypenames = NULL;
  *numimpropers      = 0;
  *impropers         = NULL;
  *impropertypes     = NULL;
  *numimpropertypes  = 0;
  *impropertypenames = NULL;
  *numcterms         = 0;
  *cterms            = NULL;
  *ctermrows         = 0;
  *ctermcols         = 0;

  *numangles = js->numangles;
  *angles = js->angles;

  *numdihedrals = js->numdihedrals;
  *dihedrals = js->dihedrals;

  *numimpropers = js->numimpropers;
  *impropers = js->impropers;

  *numcterms = js->numcterms;
  *cterms = js->cterms;
  *ctermcols = 0;
  *ctermrows = 0;

  return MOLFILE_SUCCESS;
}
#else
static int read_js_angles(void *v,
               int *numangles,    int **angles,    double **angleforces,
               int *numdihedrals, int **dihedrals, double **dihedralforces,
               int *numimpropers, int **impropers, double **improperforces,
               int *numcterms,    int **cterms,
               int *ctermcols,    int *ctermrows,  double **ctermforces) {
  jshandle *js = (jshandle *)v;

  *numangles = js->numangles;
  *angles = js->angles;
  *angleforces = NULL;

  *numdihedrals = js->numdihedrals;
  *dihedrals = js->dihedrals;
  *dihedralforces = NULL;

  *numimpropers = js->numimpropers;
  *impropers = js->impropers;
  *improperforces = NULL;

  *numcterms = js->numcterms;
  *cterms = js->cterms;
  *ctermcols = 0;
  *ctermrows = 0;
  *ctermforces = NULL;

  return MOLFILE_SUCCESS;
}
#endif

#endif


#if 1 
// XXX prototypical out-of-core trajectory analysis API
static int read_js_timestep_index_offsets(void *v, int natoms, 
                                          ptrdiff_t frameindex,
                                          int firstatom, int numatoms,
                                          fio_fd *directio_fd,
                                          ptrdiff_t *startoffset,
                                          ptrdiff_t *fileoffset,
                                          ptrdiff_t *readlen) {
  jshandle *js = (jshandle *)v;
  fio_size_t framelen;

#if JSMAJORVERSION > 1
  /* If we haven't yet read (or skipped) the structure data, then we    */
  /* need to begin by skipping past it before we try to read the        */
  /* first timestep.  In the case of files with block-aligned timesteps,*/
  /* this will also get our file pointer to the right block-aligned     */
  /* location.                                                          */
  if (!js->parsed_structure)
    read_js_structure(v, NULL, NULL);
#endif

  /* compute total read/seek size of timestep */
  framelen = js->ts_crd_padsz + js->ts_ucell_padsz;

  if (directio_fd != NULL)
    *directio_fd = js->directio_fd;

  /* compute file offset for requested timestep */
  if (fileoffset != NULL)
    *fileoffset = (frameindex * framelen) + js->ts_file_offset;

  /* compute startoffset for first requested atom */
  if (startoffset != NULL)
    *startoffset = firstatom * 3L * sizeof(float);
 
  /* compute required read size */ 
  if (readlen != NULL)
    *readlen = framelen;

  return MOLFILE_SUCCESS;
}


#if 0
static int read_js_timestep_index(void *v, int natoms, 
                                  ptrdiff_t frameindex,
                                  molfile_timestep_t *ts) {
}
#endif

#endif



static int read_js_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  jshandle *js = (jshandle *)v;
  fio_size_t framelen;

#if JSMAJORVERSION > 1
  /* If we haven't yet read (or skipped) the structure data, then we    */
  /* need to begin by skipping past it before we try to read the        */
  /* first timestep.  In the case of files with block-aligned timesteps,*/
  /* this will also get our file pointer to the right block-aligned     */
  /* location.                                                          */
  if (!js->parsed_structure)
    read_js_structure(v, NULL, NULL);
#endif

  /* compute total read/seek size of timestep */
  framelen = js->ts_crd_padsz + js->ts_ucell_padsz;

  /* if we have a valid ts pointer, read the timestep, otherwise skip it */ 
  if (ts != NULL) {
    fio_size_t readlen=0;
    fio_iovec iov[2];

    /* set unit cell pointer to the TS block-aligned buffer area */
    double *unitcell = (double *) js->directio_ucell_blkbuf;

    unitcell[0] = unitcell[2] = unitcell[5] = 1.0f;
    unitcell[1] = unitcell[3] = unitcell[4] = 90.0f;

#if defined(ENABLEJSSHORTREADS)
    /* test code for an implementation that does short reads that */
    /* skip bulk solvent, useful for faster loading of very large */
    /* structures                                                 */
    if (getenv("VMDJSMAXATOMIDX") != NULL) {
      fio_size_t bszmask;
      ptrdiff_t maxatompadsz, skipatompadsz;

      ptrdiff_t maxatomidx = atoi(getenv("VMDJSMAXATOMIDX"));
      if (maxatomidx < 0)
        maxatomidx = 0;
      if (maxatomidx >= js->natoms)
        maxatomidx = js->natoms - 1;

      /* pad max read to the start of the next block  */
      bszmask = js->directio_block_size - 1;
      maxatompadsz = ((maxatomidx*3L*sizeof(float)) + bszmask) & (~bszmask);
      skipatompadsz = js->ts_crd_padsz - maxatompadsz;

      readlen=0;
      if (js->directio_enabled) {
        if (fio_fread(ts->coords, maxatompadsz, 1, js->directio_fd) == 1)
          readlen = maxatompadsz;
        if (fio_fseek(js->directio_fd, skipatompadsz, FIO_SEEK_CUR) == 0)
          readlen += skipatompadsz;
        if (fio_fread(unitcell, js->ts_ucell_padsz, 1, js->directio_fd) == 1)
          readlen += js->ts_ucell_padsz;
      } else {
        if (fio_fread(ts->coords, maxatompadsz, 1, js->fd) == 1)
          readlen = maxatompadsz;
        if (fio_fseek(js->fd, skipatompadsz, FIO_SEEK_CUR) == 0)
          readlen += skipatompadsz;
        if (fio_fread(unitcell, js->ts_ucell_padsz, 1, js->fd) == 1)
          readlen += js->ts_ucell_padsz;
      }

#if 0
      /* clear all non-read atom coords to zeros */
      memset(ts->coords+3L*maxatomidx,0,3L*sizeof(float)*(js->natoms-maxatomidx));
#endif

    }  else {
#endif
 
    /* setup the I/O vector */
    iov[0].iov_base = (fio_caddr_t) ts->coords;   /* read coordinates    */
    iov[1].iov_base = (fio_caddr_t) unitcell;     /* read PBC unit cell  */

    if (js->directio_enabled) {
      iov[0].iov_len  = js->ts_crd_padsz;
      iov[1].iov_len  = js->ts_ucell_padsz;
    } else {
      iov[0].iov_len  = js->ts_crd_sz;
      iov[1].iov_len  = js->ts_ucell_sz;
    }
   
#if 1
    /* Use fall-back code instead of readv():                            */
    /*  Some platforms implement readv() as user level code in libc,     */
    /*  and due to POSIX atomicity requirements for readv()/writev(),    */
    /*  they may copy data to internal temp buffers, which can kill      */
    /*  performance, and in cases when doing single I/O ops on large,    */
    /*  buffers, e.g. > 2GB, can fail with shorts reads or writes...     */
    /*  On such platforms it is best to avoid using readv()/writev()...  */
    {
      int readcnt = 0;
      readlen = 0;
      if (js->directio_enabled) {
        readcnt =  fio_fread(iov[0].iov_base, iov[0].iov_len, 1, js->directio_fd);
        readcnt += fio_fread(iov[1].iov_base, iov[1].iov_len, 1, js->directio_fd);
      } else {
        fio_size_t seeklen=0;
      
        readcnt =  fio_fread(iov[0].iov_base, iov[0].iov_len, 1, js->fd);
        seeklen = js->ts_crd_padsz - js->ts_crd_sz;
        if (seeklen > 0)
          fio_fseek(js->fd, seeklen, FIO_SEEK_CUR);
        readcnt += fio_fread(iov[1].iov_base, iov[1].iov_len, 1, js->fd);
        seeklen = js->ts_ucell_padsz - js->ts_ucell_sz;
        if (seeklen > 0)
          fio_fseek(js->fd, seeklen, FIO_SEEK_CUR);
      }

      /* if both records read correctly, then the reads are okay */
      if (readcnt == 2)
        readlen = framelen;
    }
#else
    /* Do all of the reads with a single syscall, for peak efficiency.   */
    /* On smart kernels, readv() causes only one context switch, and     */
    /* can effeciently scatter the reads to the various buffers.         */
    if (js->directio_enabled) {
      readlen = fio_readv(js->directio_fd, &iov[0], 2); 
    } else {
      // XXX we can't use readv() when not using direct I/O since we 
      // can't make intervening seek calls required if the caller
      // doesn't provide appropriate buffers.
      // readlen = fio_readv(js->fd, &iov[0], 2); 

      fio_size_t seeklen=0;
      readcnt =  fio_fread(iov[0].iov_base, iov[0].iov_len, 1, js->fd);
      seeklen = js->ts_crd_padsz - js->ts_crd_sz;
      if (seeklen > 0)
        fio_fseek(js->fd, seeklen, FIO_SEEK_CUR);
      readcnt += fio_fread(iov[1].iov_base, iov[1].iov_len, 1, js->fd);
      seeklen = js->ts_ucell_padsz - js->ts_ucell_sz;
      if (seeklen > 0)
        fio_fseek(js->fd, seeklen, FIO_SEEK_CUR);
    }
#endif

#if defined(ENABLEJSSHORTREADS)
   }
#endif 
 
    /* check the number of read bytes versus what we expected */
    if (readlen != framelen) {
      if (readlen < 0) {
        perror("jsplugin) fio_readv(): ");
      } else if (readlen != 0) {
        printf("jsplugin) mismatched read: %td, expected %td\n", 
               (ptrdiff_t) readlen, (ptrdiff_t) framelen);
      }

      return MOLFILE_EOF;
    }

    /* perform byte swapping if necessary */
    if (js->reverseendian) {
      swap4_aligned(ts->coords, js->natoms * 3L);
      swap8_aligned(unitcell, 6);
    }

    /* copy unit cell values into VMD */
    ts->A = unitcell[0];
    ts->B = unitcell[1];
    ts->C = unitcell[2];
    ts->alpha = 90.0 - asin(unitcell[3]) * 90.0 / M_PI_2;
    ts->beta  = 90.0 - asin(unitcell[4]) * 90.0 / M_PI_2;
    ts->gamma = 90.0 - asin(unitcell[5]) * 90.0 / M_PI_2;
  } else {
    /* skip this frame, seek to the next frame */
    if (js->directio_enabled) {
      if (fio_fseek(js->directio_fd, framelen, FIO_SEEK_CUR)) 
        return MOLFILE_EOF;
    } else {
      if (fio_fseek(js->fd, framelen, FIO_SEEK_CUR)) 
        return MOLFILE_EOF;
    }
  }
 
  return MOLFILE_SUCCESS;
}


static void close_js_read(void *v) {
  jshandle *js = (jshandle *)v;
  fio_fclose(js->fd);

#if JSMAJORVERSION > 1
  if (js->path)
    free(js->path);

  if (js->directio_enabled)
    fio_fclose(js->directio_fd);

  if (js->directio_ucell_ptr)
    free(js->directio_ucell_ptr);

  if (js->bondfrom)
    free(js->bondfrom);
  if (js->bondto)
    free(js->bondto);
  if (js->bondorders)
    free(js->bondorders);

  /* free angle data */
  if (js->angles != NULL)
    free(js->angles);
  if (js->dihedrals != NULL)
    free(js->dihedrals);
  if (js->impropers != NULL)
    free(js->impropers);
  if (js->cterms)
    free(js->cterms);
#endif

  free(js);
}


static void *open_js_write(const char *path, const char *filetype, int natoms) {
  jshandle *js;

  js = (jshandle *) malloc(sizeof(jshandle));
  memset(js, 0, sizeof(jshandle));
#if JSMAJORVERSION > 1
  js->parsed_structure=0;
  js->directio_block_size=1;
  js->directio_ucell_ptr = NULL;
  js->directio_ucell_blkbuf = NULL;

  js->directio_enabled=0;
  js->ts_file_offset=0;
  js->ts_crd_sz=0;
  js->ts_ucell_sz=0;
  js->ts_crd_padsz=0;
  js->ts_ucell_padsz=0;
#endif

  if (fio_open(path, FIO_WRITE, &js->fd) < 0) {
    printf("jsplugin) Could not open file %s for writing\n", path);
    free(js);
    return NULL;
  }

  js->natoms = natoms;
  js->with_unitcell = 1;

  /* emit header information */
  fio_write_str(js->fd, JSHEADERSTRING);
  fio_write_int32(js->fd, JSMAGICNUMBER);
  fio_write_int32(js->fd, JSENDIANISM);
  fio_write_int32(js->fd, JSMAJORVERSION);
  fio_write_int32(js->fd, JSMINORVERSION);

  /* write number of atoms */
  fio_write_int32(js->fd, natoms);

  /* write number of frames, to be updated later */
  js->nframes = 0;
  fio_write_int32(js->fd, js->nframes);

  return js;
}


#if JSMAJORVERSION > 1

static int write_js_structure(void *mydata, int optflags,
                              const molfile_atom_t *atoms) {
  jshandle *js = (jshandle *) mydata;
  ptrdiff_t i;

  /* use block-based I/O by default when writing structures larger */
  /* than JSBLOCKIO_THRESH atoms, or when directed by the user     */
  js_blockio_check_and_set(js);

  js->optflags |= JSOPT_STRUCTURE;

  if (optflags & MOLFILE_OCCUPANCY)
    js->optflags |= JSOPT_OCCUPANCY;

  if (optflags & MOLFILE_BFACTOR)
    js->optflags |= JSOPT_BFACTOR;

  if (optflags & MOLFILE_BFACTOR)
    js->optflags |= JSOPT_BFACTOR;

  if (optflags & MOLFILE_MASS)
    js->optflags |= JSOPT_MASS;

  if (optflags & MOLFILE_CHARGE)
    js->optflags |= JSOPT_CHARGE;
 
  if (optflags & MOLFILE_RADIUS)
    js->optflags |= JSOPT_RADIUS;

  if (optflags & MOLFILE_ATOMICNUMBER)
    js->optflags |= JSOPT_ATOMICNUMBER;

  /* write flags data to the file */
  fio_write_int32(js->fd, js->optflags); 
printf("jsplugin) writing option flags: %0x08x\n", js->optflags);

  /* Check to see if block-based trajectory I/O is used  */
  /* and write out the block size for this file.         */
  if (js->optflags & JSOPT_TS_BLOCKIO) {
    fio_fwrite(&js->directio_block_size, sizeof(int), 1, js->fd);
    printf("jsplugin) Block-based I/O enabled: block size %d bytes\n", 
           js->directio_block_size);
  }

printf("jsplugin) writing structure...\n");
  /* determine whether or not this file contains structure info or not */
  if (js->optflags & JSOPT_STRUCTURE) {
    int numatomnames, numatomtypes, numresnames, numsegids, numchains;
    char **atomnames = NULL;
    char **atomtypes = NULL;
    char **resnames = NULL;
    char **segids = NULL;
    char **chains = NULL;
    short *shortbuf = NULL; /* temp buf for encoding atom records */
    int *intbuf = NULL;     /* temp buf for encoding atom records */
    float *fltbuf = NULL;   /* temp buf for encoding atom records */

    hash_t tmphash;         /* temporary hash table */
    hash_t atomnamehash;
    hash_t atomtypehash;
    hash_t resnamehash;
    hash_t segidhash;
    hash_t chainhash;
    int hashcnt;


printf("jsplugin) counting atom names, types, etc...\n");
    /* generate hash tables to count the number of unique strings */
    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].name, 0);
    numatomnames = hash_entries(&tmphash);
    hash_destroy(&tmphash);

    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].type, 0);
    numatomtypes = hash_entries(&tmphash);
    hash_destroy(&tmphash);

    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].resname, 0);
    numresnames = hash_entries(&tmphash);
    hash_destroy(&tmphash);

    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].segid, 0);
    numsegids = hash_entries(&tmphash);
    hash_destroy(&tmphash);

    hash_init(&tmphash, 127);
    for (i=0; i<js->natoms; i++)
      hash_insert(&tmphash, atoms[i].chain, 0);
    numchains = hash_entries(&tmphash);
    hash_destroy(&tmphash);
 
printf("jsplugin) writing unique string counts...\n");
printf("jsplugin) %d %d %d %d %d\n",
       numatomnames, numatomtypes, numresnames, numsegids, numchains);

    /* write block of name string table sizes */
    fio_write_int32(js->fd, numatomnames); 
    fio_write_int32(js->fd, numatomtypes); 
    fio_write_int32(js->fd, numresnames);
    fio_write_int32(js->fd, numsegids);
    fio_write_int32(js->fd, numchains); 

printf("jsplugin) writing string tables...\n");

    atomnames = (char **) malloc(numatomnames * sizeof(char *));
    atomtypes = (char **) malloc(numatomtypes * sizeof(char *));
    resnames = (char **) malloc(numresnames * sizeof(char *));
    segids = (char **) malloc(numsegids * sizeof(char *));
    chains = (char **) malloc(numchains * sizeof(char *));

printf("jsplugin)   atom names...\n");
    /* generate and write out the string tables */
    hash_init(&atomnamehash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&atomnamehash, atoms[i].name, hashcnt) == HASH_FAIL) {
        atomnames[hashcnt] = (char *) calloc(1, 16L * sizeof(char));
        strcpy(atomnames[hashcnt], atoms[i].name);
        hashcnt++;
      }
    }
    for (i=0; i<numatomnames; i++) {
      fio_fwrite(atomnames[i], 16L * sizeof(char), 1, js->fd);
    }


printf("jsplugin)   atom types...\n");
    hash_init(&atomtypehash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&atomtypehash, atoms[i].type, hashcnt) == HASH_FAIL) {
        atomtypes[hashcnt] = (char *) calloc(1, 16L * sizeof(char));
        strcpy(atomtypes[hashcnt], atoms[i].type);
        hashcnt++;
      }
    }
    for (i=0; i<numatomtypes; i++) {
      fio_fwrite(atomtypes[i], 16L * sizeof(char), 1, js->fd);
    }


printf("jsplugin)   residue names...\n");
    hash_init(&resnamehash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&resnamehash, atoms[i].resname, hashcnt) == HASH_FAIL) {
        resnames[hashcnt] = (char *) calloc(1, 8L * sizeof(char));
        strcpy(resnames[hashcnt], atoms[i].resname);
        hashcnt++;
      }
    }
    for (i=0; i<numresnames; i++) {
      fio_fwrite(resnames[i], 8L * sizeof(char), 1, js->fd);
    }


printf("jsplugin)   segment names...\n");
    hash_init(&segidhash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&segidhash, atoms[i].segid, hashcnt) == HASH_FAIL) {
        segids[hashcnt] = (char *) calloc(1, 8L * sizeof(char));
        strcpy(segids[hashcnt], atoms[i].segid);
        hashcnt++;
      }
    }
    for (i=0; i<numsegids; i++) {
      fio_fwrite(segids[i], 8L * sizeof(char), 1, js->fd);
    }


printf("jsplugin)   chain names...\n");
    hash_init(&chainhash, 127);
    for (hashcnt=0,i=0; i<js->natoms; i++) {
      /* add a new string table entry for hash inserts that don't yet exist */
      if (hash_insert(&chainhash, atoms[i].chain, hashcnt) == HASH_FAIL) {
        chains[hashcnt] = (char *) calloc(1, 2L * sizeof(char));
        strcpy(chains[hashcnt], atoms[i].chain);
        hashcnt++;
      }
    }
    for (i=0; i<numchains; i++) {
      fio_fwrite(chains[i], 2L * sizeof(char), 1, js->fd);
    }


printf("jsplugin) writing numeric field tables...\n");
    /* write out all of the atom fields */
    shortbuf = (short *) malloc(js->natoms * sizeof(short));

    /* write out atom names */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&atomnamehash, atoms[i].name);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);

    /* write out atom types */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&atomtypehash, atoms[i].type);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);

    /* write out resnames */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&resnamehash, atoms[i].resname);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);
    
    /* write out segids */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&segidhash, atoms[i].segid);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);

    /* write out chains */
    for (i=0; i<js->natoms; i++) {
      shortbuf[i] = hash_lookup(&chainhash, atoms[i].chain);
    }    
    fio_fwrite(shortbuf, js->natoms * sizeof(short), 1, js->fd);

    if (shortbuf != NULL) {
      free(shortbuf);
      shortbuf=NULL;
    }

    /* done with hash tables */
    hash_destroy(&atomnamehash);
    hash_destroy(&atomtypehash);
    hash_destroy(&resnamehash);
    hash_destroy(&segidhash);
    hash_destroy(&chainhash);


    /* 
     * write out integer data blocks 
     */
    intbuf = (int *) malloc(js->natoms * sizeof(int));

printf("jsplugin)   residue indices...\n");
    /* write out resid */
    for (i=0; i<js->natoms; i++) {
      intbuf[i] = atoms[i].resid;
    }    
    fio_fwrite(intbuf, js->natoms * sizeof(int), 1, js->fd);
     
    if (intbuf != NULL) {
      free(intbuf);
      intbuf = NULL;
    }

printf("jsplugin) writing optional per-atom tables...\n");
    /*
     * write out optional single-precision float data blocks
     */ 
    if (js->optflags & (JSOPT_OCCUPANCY | JSOPT_BFACTOR | 
        JSOPT_MASS | JSOPT_RADIUS | JSOPT_CHARGE)) 
      fltbuf = (float *) malloc(js->natoms * sizeof(float));

    /* write out optional data if it exists */

    if (js->optflags & JSOPT_OCCUPANCY) {
printf("jsplugin)   writing occupancy...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].occupancy;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (js->optflags & JSOPT_BFACTOR) {
printf("jsplugin)   writing bfactor...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].bfactor;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (js->optflags & JSOPT_MASS) { 
printf("jsplugin)   writing mass...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].mass;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (js->optflags & JSOPT_CHARGE) { 
printf("jsplugin)   writing charge...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].charge;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (js->optflags & JSOPT_RADIUS) { 
printf("jsplugin)   writing radius...\n");
      for (i=0; i<js->natoms; i++) {
        fltbuf[i] = atoms[i].radius;
      }    
      fio_fwrite(fltbuf, js->natoms * sizeof(float), 1, js->fd);
    }

    if (fltbuf != NULL) {
      free(fltbuf);
      fltbuf=NULL;
    }


    /*
     * write out optional integer data blocks
     */ 
    if (js->optflags & JSOPT_ATOMICNUMBER)
      intbuf = (int *) malloc(js->natoms * sizeof(int));

    if (js->optflags & JSOPT_ATOMICNUMBER) { 
printf("jsplugin)   writing atomic number...\n");
      for (i=0; i<js->natoms; i++) {
        intbuf[i] = atoms[i].atomicnumber;
      }    
      fio_fwrite(intbuf, js->natoms * sizeof(int), 1, js->fd);
    }

    if (intbuf != NULL) {
      free(intbuf);
      intbuf = NULL;
    }


    /*
     * write out bonds and fractional bond orders
     */ 
    if (js->optflags & JSOPT_BONDS) {
printf("jsplugin) writing bonds...\n");
      fio_fwrite(&js->nbonds, sizeof(int), 1, js->fd);
      fio_fwrite(js->bondfrom, js->nbonds * sizeof(int), 1, js->fd);
      fio_fwrite(js->bondto, js->nbonds * sizeof(int), 1, js->fd);

      if (js->optflags & JSOPT_BONDORDERS) {
printf("jsplugin) writing bond orders...\n");
        fio_fwrite(js->bondorders, js->nbonds * sizeof(float), 1, js->fd);
      }
    }

    /*
     * write out angles/dihedrals/impropers/cross-terms
     */
    if (js->optflags & JSOPT_ANGLES) {
printf("jsplugin) writing angles/dihedrals/impropers...\n");
      fio_fwrite(&js->numangles, sizeof(int), 1, js->fd);
      fio_fwrite(js->angles, sizeof(int)*3L*js->numangles, 1, js->fd);

      fio_fwrite(&js->numdihedrals, sizeof(int), 1, js->fd);
      fio_fwrite(js->dihedrals, sizeof(int)*4L*js->numdihedrals, 1, js->fd);

      fio_fwrite(&js->numimpropers, sizeof(int), 1, js->fd);
      fio_fwrite(js->impropers, sizeof(int)*4L*js->numimpropers, 1, js->fd);
    }
    if (js->optflags & JSOPT_CTERMS) {
printf("jsplugin) writing cross-terms\n");
      fio_fwrite(&js->numcterms, sizeof(int), 1, js->fd);
      fio_fwrite(js->cterms, sizeof(int)*8L*js->numcterms, 1, js->fd);
    }

    /* update the file offset for the first timestep */
    js_calc_timestep_blocking_info(js);

    return MOLFILE_SUCCESS;
  }

  /* update the file offset for the first timestep */
  js_calc_timestep_blocking_info(js);

  /* else, we have no structure information */
  return MOLFILE_NOSTRUCTUREDATA;
}


static int write_js_bonds(void *mydata, int nbonds, int *fromptr, int *toptr, 
                          float *bondorder,  int *bondtype, 
                          int nbondtypes, char **bondtypename) {
  jshandle *js = (jshandle *) mydata;

#if defined(INFOMSGS)
    if (js->verbose) {
      printf("jsplugin) write_js_bonds():\n");
      printf("jsplugin) storing bond info for writing...\n");
      printf("jsplugin) %d %d\n", nbonds, nbondtypes);
    }
#endif

  if (nbonds > 0 && fromptr != NULL && toptr != NULL) {
    js->optflags |= JSOPT_BONDS; 

    /* save bond info until we actually write out the structure file */
    js->nbonds = nbonds;
    js->bondfrom = (int *) malloc(nbonds * sizeof(int));
    memcpy(js->bondfrom, fromptr, nbonds * sizeof(int));
    js->bondto = (int *) malloc(nbonds * sizeof(int));
    memcpy(js->bondto, toptr, nbonds * sizeof(int));

    if (bondorder != NULL) {
      js->optflags |= JSOPT_BONDORDERS;
      js->bondorders = (float *) malloc(nbonds * sizeof(float));
      memcpy(js->bondorders, bondorder, nbonds * sizeof(float));
    }
  }

  return MOLFILE_SUCCESS;
}

#if vmdplugin_ABIVERSION > 14
static int write_js_angles(void * v, int numangles, const int *angles,
                           const int *angletypes, int numangletypes,
                           const char **angletypenames, int numdihedrals, 
                           const int *dihedrals, const int *dihedraltype,
                           int numdihedraltypes, const char **dihedraltypenames,
                           int numimpropers, const int *impropers, 
                           const int *impropertypes, int numimpropertypes, 
                           const char **impropertypenames, int numcterms, 
                           const int *cterms, int ctermcols, int ctermrows) {
  jshandle *js = (jshandle *) v;

  /* save info until we actually write out the structure file */
  js->numangles = numangles;
  js->numdihedrals = numdihedrals;
  js->numimpropers = numimpropers;
  js->numcterms = numcterms;

#if defined(INFOMSGS)
  if (js->verbose) {
    printf("jsplugin) write_js_angles():\n");
    printf("jsplugin) storing angles/dihedrals/impropers for writing...\n");
    printf("jsplugin) %d %d %d %d\n",
           numangles, numdihedrals, numimpropers, numcterms);
  }
#endif

  if (js->numangles > 0 || js->numdihedrals > 0 || js->numimpropers > 0) {
    js->optflags |= JSOPT_ANGLES;

    js->angles = (int *) malloc(3L*js->numangles*sizeof(int));
    memcpy(js->angles, angles, 3L*js->numangles*sizeof(int));
    js->dihedrals = (int *) malloc(4L*js->numdihedrals*sizeof(int));
    memcpy(js->dihedrals, dihedrals, 4L*js->numdihedrals*sizeof(int));
    js->impropers = (int *) malloc(4L*js->numimpropers*sizeof(int));
    memcpy(js->impropers, impropers, 4L*js->numimpropers*sizeof(int));
  }
  if (js->numcterms > 0) {
    js->optflags |= JSOPT_CTERMS;

    js->cterms = (int *) malloc(8L*js->numcterms*sizeof(int));
    memcpy(js->cterms, cterms, 8L*js->numcterms*sizeof(int));
  }

  return MOLFILE_SUCCESS;
}
#else
static int write_js_angles(void * v,
        int numangles,    const int *angles,    const double *angleforces,
        int numdihedrals, const int *dihedrals, const double *dihedralforces,
        int numimpropers, const int *impropers, const double *improperforces,
        int numcterms,   const int *cterms,
        int ctermcols, int ctermrows, const double *ctermforces) {
  jshandle *js = (jshandle *) v;

  /* save info until we actually write out the structure file */
  js->numangles = numangles;
  js->numdihedrals = numdihedrals;
  js->numimpropers = numimpropers;
  js->numcterms = numcterms;

  if (js->numangles > 0 || js->numdihedrals > 0 || js->numimpropers > 0) {
    js->optflags |= JSOPT_ANGLES;

    js->angles = (int *) malloc(3L*js->numangles*sizeof(int));
    memcpy(js->angles, angles, 3L*js->numangles*sizeof(int));
    js->dihedrals = (int *) malloc(4L*js->numdihedrals*sizeof(int));
    memcpy(js->dihedrals, dihedrals, 4L*js->numdihedrals*sizeof(int));
    js->impropers = (int *) malloc(4L*js->numimpropers*sizeof(int));
    memcpy(js->impropers, impropers, 4L*js->numimpropers*sizeof(int));
  }
  if (js->numcterms > 0) {
    js->optflags |= JSOPT_CTERMS;

    js->cterms = (int *) malloc(8L*js->numcterms*sizeof(int));
    memcpy(js->cterms, cterms, 8L*js->numcterms*sizeof(int));
  }

  return MOLFILE_SUCCESS;
}
#endif
#endif


static int write_js_timestep(void *v, const molfile_timestep_t *ts) { 
  jshandle *js = (jshandle *)v;
  double *unitcell=NULL;
  ptrdiff_t zeropadsz=0;

  /* If no structure data was written and this is the first timestep */
  /* we must complete writing the file header and performing the     */
  /* seek to the next filesystem block and VM-page boundary when     */
  /* using direct I/O APIs...                                        */
  if (js->directio_ucell_blkbuf == NULL) {
    printf("jsplugin) no structure data, writing timesteps only...\n");

    /* use block-based I/O by default when writing structures larger */
    /* than JSBLOCKIO_THRESH atoms, or when directed by the user     */
    js_blockio_check_and_set(js);

    /* write flags data to the file */
    fio_write_int32(js->fd, js->optflags); 
    printf("jsplugin) writing option flags: %0x08x\n", js->optflags);

    /* Check to see if block-based trajectory I/O is used  */
    /* and write out the block size for this file.         */
    if (js->optflags & JSOPT_TS_BLOCKIO) {
      fio_fwrite(&js->directio_block_size, sizeof(int), 1, js->fd);
      printf("jsplugin) Block-based I/O enabled: block size %d bytes\n", 
             js->directio_block_size);
    }

    /* update the file offset for the first timestep */
    js_calc_timestep_blocking_info(js);
  }

  /* set unit cell pointer to the TS block-aligned buffer area */
  unitcell = (double *) js->directio_ucell_blkbuf;

  js->nframes++; /* increment frame count written to the file so far */

  unitcell[0] = ts->A;
  unitcell[1] = ts->B;
  unitcell[2] = ts->C;
  unitcell[3] = sin((M_PI_2 / 90.0) * (90.0 - ts->alpha));
  unitcell[4] = sin((M_PI_2 / 90.0) * (90.0 - ts->beta));
  unitcell[5] = sin((M_PI_2 / 90.0) * (90.0 - ts->gamma));

  /* coordinates for all atoms */
  if (fio_fwrite(ts->coords, js->ts_crd_sz, 1, js->fd) != 1) {
    printf("jsplugin) Error writing timestep coords!\n");
    return MOLFILE_ERROR;
  }

  /* correctly handle block-based direct-I/O output format       */
  /* write out coord padding bytes using zero buffer in jshandle */
  zeropadsz = js->ts_crd_padsz - js->ts_crd_sz;
  if (zeropadsz > 0) {
    if ((zeropadsz > MOLFILE_DIRECTIO_MAX_BLOCK_SIZE) ||
        (fio_fwrite(js->blockpad, zeropadsz, 1, js->fd) != 1)) {
      printf("jsplugin) Error writing timestep coord padding!\n");
      return MOLFILE_ERROR;
    } 
  }

  /* PBC unit cell info */ 
  if (fio_fwrite(unitcell, js->ts_ucell_sz, 1, js->fd) != 1) {
    printf("jsplugin) Error writing timestep unit cell!\n");
    return MOLFILE_ERROR;
  }

  /* correctly handle block-based direct-I/O output format       */
  /* write out PBC padding bytes using zero buffer in jshandle */
  zeropadsz = js->ts_ucell_padsz - js->ts_ucell_sz;
  if (zeropadsz > 0) {
    if ((zeropadsz > MOLFILE_DIRECTIO_MAX_BLOCK_SIZE) ||
        (fio_fwrite(js->blockpad, zeropadsz, 1, js->fd) != 1)) {
      printf("jsplugin) Error writing timestep PBC padding!\n");
      return MOLFILE_ERROR;
    } 
  }

  return MOLFILE_SUCCESS;
}


static void close_js_write(void *v) {
  jshandle *js = (jshandle *)v;

  /* update the trajectory header information */
  fio_fseek(js->fd, JSNFRAMESOFFSET, FIO_SEEK_SET);
  fio_write_int32(js->fd, js->nframes);
  fio_fseek(js->fd, 0, FIO_SEEK_END);

  fio_fclose(js->fd);

#if JSMAJORVERSION > 1
  if (js->directio_ucell_ptr)
    free(js->directio_ucell_ptr);

  if (js->bondfrom)
    free(js->bondfrom);
  if (js->bondto)
    free(js->bondto);
  if (js->bondorders)
    free(js->bondorders);

  if (js->angles)
    free(js->angles);
  if (js->dihedrals)
    free(js->dihedrals);
  if (js->impropers)
    free(js->impropers);
  if (js->cterms)
    free(js->cterms);
#endif

  free(js);
}


/*
 * Initialization stuff here
 */
static molfile_plugin_t plugin;

#if !defined(VMDJSPLUGININCLUDESRC)

VMDPLUGIN_API int VMDPLUGIN_init() {
  memset(&plugin, 0, sizeof(molfile_plugin_t));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "js";
  plugin.prettyname = "js";
  plugin.author = "John Stone";
  plugin.majorv = JSMAJORVERSION;
  plugin.minorv = JSMINORVERSION;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
  plugin.filename_extension = "js";
  plugin.open_file_read = open_js_read;
#if JSMAJORVERSION > 1
  plugin.read_structure = read_js_structure;
  plugin.read_bonds = read_js_bonds;
  plugin.read_angles = read_js_angles;
#endif
  plugin.read_next_timestep = read_js_timestep;
  plugin.close_file_read = close_js_read;
  plugin.open_file_write = open_js_write;
#if JSMAJORVERSION > 1
  plugin.write_structure = write_js_structure;
  plugin.write_bonds = write_js_bonds;
  plugin.write_angles = write_js_angles;
#endif
  plugin.write_timestep = write_js_timestep;
  plugin.close_file_write = close_js_write;
#if vmdplugin_ABIVERSION > 17
  plugin.read_timestep_pagealign_size = read_js_timestep_pagealign_size;
#endif
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb) {
  (*cb)(v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini() {
  return VMDPLUGIN_SUCCESS;
}

#endif
  
#ifdef TEST_JSPLUGIN

#include <sys/time.h>

#if defined(ENABLECUDATESTS)
#include <cuda_runtime.h>

#if defined(ENABLECUDAGDS)
#include <cufile.h>
#endif
#endif

/* get the time of day from the system clock, and store it (in seconds) */
double time_of_day(void) {
#if defined(_MSC_VER)
  double t;

  t = GetTickCount();
  t = t / 1000.0;

  return t;
#else
  struct timeval tm;
  struct timezone tz;

  gettimeofday(&tm, &tz);
  return((double)(tm.tv_sec) + (double)(tm.tv_usec)/1000000.0);
#endif
}

int main(int argc, char *argv[]) {
  molfile_timestep_t timestep;
  float *coords0=NULL, *aligncoords0=NULL;
  float *coords1=NULL, *aligncoords1=NULL;

  void *v;
  jshandle *js;
  int natoms, i;
  ptrdiff_t sz, blocksz;
  float sizeMB =0.0, totalMB = 0.0;
  double starttime, endtime, totaltime = 0.0;
  int do_io = 1;
  int verbose = 0;
  int overlapiogpu = 1;

  printf("Standalone tests for JS plugin:\n");
  
  if (getenv("VMDJSNOIO") != NULL)
    do_io = 0;

  if (getenv("VMDJSVERBOSE") != NULL || getenv("VMDJSTESTVERBOSE"))
    verbose = 1;    

  if (do_io)
    printf("  Timestep disk I/O enabled.\n");
  else
    printf("  Timestep disk I/O DISABLED.\n");

#if defined(ENABLECUDATESTS)
  printf("  CUDA GPU support compiled in.\n");

  // If the code is compiled with CUDA support, we benchmark 
  // host I/O immediately followed by host-GPU copies of each timestep
  cudaError_t crc;
  cudaStream_t devstream;
  ptrdiff_t maxatomidx=-1;
  int devcount;
  float *devptr=NULL;

  crc = cudaGetDeviceCount(&devcount);
  printf("  GPU device count: %d\n", devcount);
  if (devcount==0)
    printf("  No GPU devices, continuing with host only...\n");

  // Only do the CUDA tests if asked to
  if (getenv("VMDJSCUDATESTS") == NULL) {
    devcount = 0;
    printf("  GPU tests disabled.\n");
    printf("  Enable GPU tests with VMDJSCUDATESTS env variable\n");
  } else {
    printf("  Disable GPU tests by unsetting VMDJSCUDATESTS env variable\n");
  }

#if defined(ENABLEJSSHORTREADS)
  /* test code for an implementation that does short reads that */
  /* skip bulk solvent, useful for faster loading of very large */
  /* structures                                                 */
  if (getenv("VMDJSMAXATOMIDX") != NULL) {
    fio_size_t bszmask;

    maxatomidx = atoi(getenv("VMDJSMAXATOMIDX"));
    if (maxatomidx < 0)
      maxatomidx = 0;
    if (maxatomidx >= js->natoms)
      maxatomidx = js->natoms - 1;

    printf("jsplugin) Short-copies of GPU timesteps enabled: %ld / %ld atoms (%.2f%%)\n",
           maxatomidx, js->natoms, 100.0*(maxatomidx+1) / ((float) js->natoms));
  }
#endif
#endif


  while (--argc) {
    int syncframe;
    ++argv; 
    natoms = 0;
    v = open_js_read(*argv, "js", &natoms);
    if (!v) {
      printf("jsplugin) open_js_read failed for file %s\n", *argv);
      return 1;
    }
    js = (jshandle *)v;
    sizeMB = ((natoms * 3.0) * js->nframes * 4.0) / (1024.0 * 1024.0);
    totalMB += sizeMB; 
    printf("jsplugin) file: %s\n", *argv);
    printf("jsplugin)   %d atoms, %d frames, size: %6.1fMB\n", natoms, js->nframes, sizeMB);

    starttime = time_of_day();

    /* ensure we have a large enough allocation so we can align */
    /* the starting pointer to a blocksz page boundary          */
    blocksz = MOLFILE_DIRECTIO_MIN_BLOCK_SIZE;
    sz = 3L*sizeof(float)*natoms + blocksz;

    /* pad the allocation to an even multiple of the block size */
    size_t blockpadsz = (sz + (blocksz - 1)) & (~(blocksz - 1));

    /* allocate multi buffers so we can have concurrent work/transfers/reads */
    aligncoords0 = (float *) alloc_aligned_ptr(sz, blocksz, (void**) &coords0);
    aligncoords1 = (float *) alloc_aligned_ptr(sz, blocksz, (void**) &coords1);

#if vmdplugin_ABIVERSION > 17
    int filepgalignsz = 1;
    read_js_timestep_pagealign_size(v, &filepgalignsz);
    if (filepgalignsz != blocksz) {
      printf("jsplugin) Plugin-returned page alignment size doesn't match!\n");
    } else {
      printf("jsplugin) Page alignment size: %d\n", filepgalignsz);
    }
#endif

#if defined(ENABLECUDATESTS)
#if defined(ENABLECUDAGDS)
    int cufileinuse=0;
    CUfileHandle_t cfh;
    CUfileDescr_t cfhdesc;
    CUfileError_t cferr;
    memset(&cfh, 0, sizeof(cfh));
    memset(&cfhdesc, 0, sizeof(cfhdesc));
#endif

    if (crc == cudaSuccess && devcount > 0) {
      cudaSetDevice(0);
      printf("jsplugin) allocating GPU memory buffer for CUDA tests...\n");
      crc = cudaMalloc((void**) &devptr, blockpadsz);
      if (crc != cudaSuccess) {
        printf("Failed to allocate GPU buffer!\n");
        return -1;
      }

      cudaStreamCreate(&devstream);

#if defined(ENABLECUDAGDS)
      cuFileBufRegister(devptr, blockpadsz, 0); 

      cfhdesc.handle.fd = js->directio_fd; // typedef of Unix FD
      cfhdesc.type = CU_FILE_EXTERNAL_MEMORY_HANDLE_TYPE_OPAQUE_FD;
      cferr = cuFileImportExternalFile(&cfh, &cfhdesc);
      if (cferr.err != CU_FILE_SUCCESS) {
        printf("Failed to import file handle for use by cuFile APIs!\n");
        return -1;
      }
      cufileinuse=1;
#endif
    }
#endif

    /* loop over all timesteps ... */
    for (syncframe=0,i=0; i<js->nframes; i++) {
      if (do_io) {
        /* read even/odd frame into alternating buffers so     */
        /* that we can overlap a read with an ongoing GPU copy */
        if (i & 1) 
          timestep.coords = aligncoords1;
        else  
          timestep.coords = aligncoords0;

        /* disk I/O is synchronous */
        if (verbose) {
          printf("%sreading frame[%d]...", (i!=0) ? "\r" : "", i);
          fflush(stdout);
        }
        int rc=0;
#if defined(ENABLECUDAGDS)
        if (cufileinuse) {
#if 0
          /* read an even multiple of the block size */
          ptrdiff_t rsz = natoms * 3L * sizeof(float);
          rsz = (rsz + (blocksz - 1)) & (~(blocksz - 1));
          ptrdiff_t foffset = i * rsz;
#endif
          ptrdiff_t startoffset, foffset, readlen;
          fio_fd directio_fd;
          read_js_timestep_index_offsets(v, natoms, i, 0, natoms,
                                         &directio_fd,
                                         &startoffset,
                                         &foffset,
                                         &readlen);

printf("cuFileRead(): offset %ld  readlen: %ld\n", foffset, readlen);
          ptrdiff_t ret = 0;
          ret = cuFileRead(cfh, (char *) devptr, readlen, foffset);
          if (ret < 0) {
            const char *descp = "unknown error code";
#if 0
            // XXX this requires linkage with libcuda.so, which is
            //     against convention, so we avoid it for now
            if (cuGetErrorName(ret, &descp) != CUDA_SUCCESS)
              descp = "unknown cuda error";
#endif

            printf("Error: cuFileRead(): %ld, '%s'\n", ret, descp);
            return -1;
          }
        } else
#else
          rc = read_js_timestep(v, natoms, &timestep);
#endif

        if (rc) {
          printf("jsplugin) error in read_js_timestep on frame %d\n", i);
          /* return 1; */
        }
      }

#if defined(ENABLECUDATESTS)
      if (crc == cudaSuccess && devcount > 0) {
#if defined(ENABLECUDAGDS)
        if (!cufileinuse) {
#endif
          /* allow overlap of async memcpy with next round of direct I/O */
          if (overlapiogpu) {
            if (verbose) {
              printf("sync frame[%d]...", syncframe);
              fflush(stdout);
            }
            cudaStreamSynchronize(devstream);
          }

          if (verbose) {
            printf("cudaMemcpyAsync() frame[%d]...", i);
            fflush(stdout);
          }

          size_t bsz = (maxatomidx >= 0) ? (maxatomidx+1) : natoms;
          bsz *= 3L*sizeof(float);
          crc = cudaMemcpyAsync(devptr, timestep.coords, bsz, cudaMemcpyHostToDevice, devstream);
          syncframe=i;

          if (!overlapiogpu) {
            if (verbose) {
              printf("sync frame[%d]...", syncframe);
              fflush(stdout);
            }
            cudaStreamSynchronize(devstream);
          }
  
          if (verbose) {
            printf("       ");
            fflush(stdout);
          }
#if defined(ENABLECUDAGDS)
        }
#endif
      }
#endif
    }

#if defined(ENABLECUDATESTS)
    /* wait for last GPU memcpy */
    cudaStreamSynchronize(devstream);
    printf("\n");

    /* wait for any pending GPU calls to complete */
    cudaDeviceSynchronize();

#if defined(ENABLECUDAGDS)
    if (cufileinuse) {
      cuFileBufDeregister(devptr); 
      cuFileDestroyFile(&cfh);
    }
#endif
#endif

    endtime = time_of_day();
    close_js_read(v);

#if defined(ENABLECUDATESTS)
    cudaStreamDestroy(devstream);
    if (crc == cudaSuccess && devcount > 0) {
      cudaFree(devptr);
    }
#endif
    free(coords0);
    free(coords1);

    totaltime += endtime - starttime;
    printf("jsplugin)  Time: %5.1f seconds\n", endtime - starttime);
    printf("jsplugin)  Speed: %5.1f MB/sec, %5.1f timesteps/sec\n", sizeMB / (endtime - starttime), (js->nframes / (endtime - starttime)));
  }
  printf("jsplugin) Overall Size: %6.1f MB\n", totalMB);
  printf("jsplugin) Overall Time: %6.1f seconds\n", totaltime);
  printf("jsplugin) Overall Speed: %5.1f MB/sec\n", totalMB / totaltime);
  return 0;
}
      
#endif

