#ifndef __config_h
#define __config_h

#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#ifdef DEBUG
   #define __DEBUG (DEBUG)
#else
   #define __DEBUG (1)
#endif
#ifdef PROFILE
   #define __PROFILE (PROFILE)
#else
   #define __PROFILE (1)
#endif


// Enable/Disable the assert macro ...
#ifndef DEBUG
#define NDEBUG 1
#warning 'disabled macro assert.h'
#endif

#if defined(DEBUG)
#if (DEBUG > 0)
#undef NDEBUG
//#else
//#warning 'enabled macro assert.h'
#endif
#endif

#include <assert.h>

#endif
