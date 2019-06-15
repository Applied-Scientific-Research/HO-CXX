#ifndef __has_parallel_scoped_lock_h
#define __has_parallel_scoped_lock_h

#include <utils/timer.h>
#include <splib/attributes.h>
#include <splib/spblas.h>

#include <parallel/parallel.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace parallel
{

#ifdef _OPENMP
struct MutexType
{
   MutexType() { omp_init_lock(&this->lock); }
   ~MutexType() { omp_destroy_lock(&this->lock); }
   void Lock() { omp_set_lock(&this->lock); }
   void Unlock() { omp_unset_lock(&this->lock); }
   MutexType(const MutexType& ) { omp_init_lock(&this->lock); }
   MutexType& operator= (const MutexType& ) { return *this; }
 public:
   omp_lock_t lock;
};
#else
struct MutexType
{
   void Lock() {};
   void Unlock() {};
};
#endif

struct ScopedLock
{
   explicit ScopedLock(MutexType& m) : m_mutex(m), m_locked(true) { this->m_mutex.Lock(); }
   ~ScopedLock() { this->Unlock(); }
   void Unlock() { if(not(this->m_locked)) return; this->m_locked=false; this->m_mutex.Unlock(); }
   void Relock() { if(this->m_locked) return; this->m_mutex.Lock(); this->m_locked=true; }
 private:
   MutexType& m_mutex;
   bool m_locked;
 private: // prevent copying the scoped lock.
   void operator=(const ScopedLock&);
   ScopedLock(const ScopedLock&);
};

// forward declaration ...
//class ThreadBarrier
//{
//};

} // namespace parallel

#endif
