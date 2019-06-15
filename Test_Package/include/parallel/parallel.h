#ifndef __has_parallel_parallel_h
#define __has_parallel_parallel_h

#include <utils/xml.h>

namespace parallel
{
   //enum { max_threads = 64 };
   enum { max_threads = 240 };
   //enum { max_partitions = 16 };

   //static int num_threads = 1;
   extern int num_threads;

   extern int verbosity;

   //static bool is_initialized = 0;
   extern bool is_initialized;

#ifdef ENABLE_XML2
   //static xmlDoc *xml_doc = NULL;
   extern xmlDoc *xml_doc;
#endif

   struct library_variables_t
   {
      library_variables_t (void);
      ~library_variables_t(void);
   };

   //static library_variables_t library_variables;

   void initialize(void);
   int get_num_threads(void);
   int get_thread_num(void);

} // parallel

#endif
