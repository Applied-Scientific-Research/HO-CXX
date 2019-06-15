#ifndef __has_partitioned_objects_h
#define __has_partitioned_objects_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include <utils/config.h>
#include <parallel/numa.h>
#include <parallel/templates/partitioner.inl>

#include <splib/spblas.h>
#include <splib/copy.h>
#include <splib/multiply.h>

#include <iostream>
#include <typeinfo>

namespace parallel
{

enum { EnableNestedParallelism = 0 };

#ifndef __DEFAULT_NUM_PARTITIONS
#define __DEFAULT_NUM_PARTITIONS (2)
#endif
enum { DefaultNumberPartitions = __DEFAULT_NUM_PARTITIONS };

//static bool enable_partitioning = false;
enum { enable_partitioning = 0 };

template <typename VectorType>
struct partitioned_vector
{
   //typedef typename VectorType::format_type	format_type;
   typedef splib::partitioned_vector_format	format_type;
   typedef typename VectorType::value_type	value_type;
   typedef typename VectorType::index_type	index_type;

#ifdef USE_LIBNUMA
   typedef parallel::numa_local_allocator<value_type>
						allocator_type;
#else
   typedef typename VectorType::allocator_type	allocator_type;
#endif

   typedef typename VectorType::template rebind<value_type,allocator_type>::other
						vector_type;

   int m_size;
   int num_partitions;
   vector_type *partition;
   splib::Vector<int> offset;

   //partitioned_vector (const VectorType &v, const int n = parallel::numa_params.num_nodes)
   partitioned_vector (const VectorType &v, const int n = DefaultNumberPartitions)
      :
         m_size(v.size()), num_partitions(n), partition(NULL), offset(n+1)
   {
      if (__DEBUG) std::cout << "... constructing (copy) partitioned vector ... " << typeid(*this).name() << " " << num_partitions << std::endl;

      this->alloc();
      this->copy(v);
   }

   //partitioned_vector (const int _size, const int n = parallel::numa_params.num_nodes)
   partitioned_vector (const int _size, const int n = DefaultNumberPartitions)
      :
         m_size(_size), num_partitions(n), partition(NULL), offset(n+1)
   {
      if (__DEBUG) std::cout << "... constructing (empty) partitioned vector ... " << typeid(*this).name() << " " << num_partitions << std::endl;

      this->alloc();
   }

   partitioned_vector (void)
      :
         m_size(0), num_partitions(0), partition(NULL), offset()
   {
      if (__DEBUG) std::cout << "... constructing (nil) partitioned vector ... " << typeid(*this).name() << std::endl;
   }

   void alloc (void)
   {
      if (__DEBUG) std::cout << "... allocing partitioned vector ... " << typeid(*this).name() << " " << num_partitions << std::endl;

      if (this->partition)
      {
         fprintf(stderr,"partitions already allocated %s %d\n", __FILE__, __LINE__);
         exit(-1);
      }

      partition = new vector_type [num_partitions];
      if (partition == NULL)
         throw std::bad_alloc();

      this->offset.resize(num_partitions+1);

      parallel::partitioner (&this->offset[0], this->num_partitions, 0, this->size());

#ifdef USE_LIBNUMA
      // Get current mask ...
      nodemask_t default_mask;
      if (parallel::numa_params.is_initialized)
         default_mask = numa_get_membind();
      //printf("default_mask  equal(all_nodes)=%d\n", nodemask_equal(&default_mask,&numa_all_nodes));

      if (this->num_partitions != numa_params.num_nodes)
      {
         fprintf(stderr,"this->num_partitions != numa_params.num_nodes %s %d\n", __FILE__, __LINE__);
         exit(-1);
      }
#endif

      for (int p = 0; p < this->num_partitions; ++p)
      {
#ifdef USE_LIBNUMA
         if (parallel::numa_params.is_initialized)
         {
            nodemask_t mask;
            nodemask_zero(&mask);
            nodemask_set(&mask,p);
            numa_set_membind(&mask);

            nodemask_t nodemask = numa_get_membind();
            if (__DEBUG)
            {
               printf("node_mask = %d ", p);
               for (int i = 0; i < parallel::numa_params.num_nodes; ++i)
                  printf("%1d", nodemask_isset(&nodemask,i));
               printf("\n");
            }
         }
#endif
         int start = this->offset[p];
         int stop  = this->offset[p+1];
         int nvals = stop - start;

         // Need to specify local allocation somehow ... just set local alloc policy ...
         this->partition[p].resize(nvals);

         if (__DEBUG) printf("partition[%d]: start=%d, stop=%d, nvals=%d\n", p, start, stop, this->partition[p].size());
      }

#ifdef USE_LIBNUMA
      if (parallel::numa_params.is_initialized)
         numa_set_membind(&default_mask);
#endif
   }

   ~partitioned_vector()
   {
      this->offset.clear();
      if (partition)
         delete [] partition;
   }

   int size(void) const { return m_size; }

   //void resize (const int _size, const int n = parallel::numa_params.num_nodes)
   void resize (const int _size, const int n = DefaultNumberPartitions)
   {
      this->~partitioned_vector();

      this->m_size = _size;
      this->num_partitions = n;

      if (__DEBUG) std::cout << "... resizing partitioned vector ... " << typeid(*this).name() << " " << num_partitions << std::endl;

      this->alloc();
   }

   void copy (const VectorType &v)
   {
      for (int p = 0; p < this->num_partitions; ++p)
      {
         int start = this->offset[p];
         int stop  = this->offset[p+1];
         int nvals = stop - start;

         const VectorType vtmp(&v[0] + start, nvals);

         //splib::copy(splib::make_vector_view(v, nvals, start), this->partition[p]);
         splib::spblas::copy(vtmp, this->partition[p]);
//       for (int i = 0; i < nvals; ++i)
//          this->partition[p][i] = v[start+i];
      }
   }

};

template <typename MatrixType>
struct partitioned_matrix
   : public splib::BaseMatrix
{
   enum { FormatTag = splib::FormatTags::SPARSE };
   enum { ValueTag = splib::ValueTags::TYPE< typename MatrixType::value_type>::tag };

   typedef typename MatrixType::value_type	value_type;
   typedef typename MatrixType::index_type	index_type;
#ifdef USE_LIBNUMA
   typedef parallel::numa_local_allocator<value_type> allocator_type;
#else
   typedef typename MatrixType::allocator_type	allocator_type;
#endif
   typedef splib::partitioned_matrix_format	format_type;

   typedef typename MatrixType::template rebind<value_type,allocator_type>::other
						matrix_type;

   typedef splib::BaseMatrix Parent;

   // Just an array of base matrix
   int num_partitions;
   matrix_type *partition;
   splib::Vector<int> offset;

   void _copy(const MatrixType &M, splib::csr_format)
   {
#ifdef USE_LIBNUMA
      nodemask_t default_mask;
      // Get current mask ...
      if (parallel::numa_params.is_initialized)
      {
         default_mask = numa_get_membind();
         if (__DEBUG)
         {
            printf("default_mask  equal(all_nodes)=%d\n", nodemask_equal(&default_mask,&numa_all_nodes));
            for (int i = 0; i <= numa_max_node(); ++i)
               printf("node %d on/off = %1d\n", i, nodemask_isset(&default_mask,i));
         }
      }

      if (this->num_partitions != numa_params.num_nodes)
      {
         fprintf(stderr,"this->num_partitions != numa_params.num_nodes %d %d %s %d\n", this->num_partitions, numa_params.num_nodes, __FILE__, __LINE__);
         exit(-1);
      }
#endif

      //int nops(0);

      for (int p = 0; p < this->num_partitions; ++p)
      {
#ifdef USE_LIBNUMA
         if (parallel::numa_params.is_initialized)
         {
            nodemask_t mask;
            nodemask_zero(&mask);
            nodemask_set(&mask,p);
            numa_set_membind(&mask);

            nodemask_t nodemask = numa_get_membind();
            if (__DEBUG)
            {
               printf("node_mask = %d ", p);
               for (int i = 0; i <= numa_max_node(); ++i)
                  printf("%1d", nodemask_isset(&nodemask,i));
               printf("\n");
            }
         }
#endif
         int row_start = this->offset[p];
         int row_stop  = this->offset[p+1];
         int num_rows  = row_stop - row_start;
         int num_values= M.rowptr[row_stop] - M.rowptr[row_start];

         // Need to specify local allocation somehow ... just set local alloc policy ...
         this->partition[p].resize (num_rows, M.num_columns, num_values);

         this->partition[p].rowptr[0] = 0;
         for (int i = 0; i < num_rows; ++i)
         {
            index_type global_k = M.rowptr[row_start+i];
            index_type nvals = M.rowptr[row_start+i+1] - M.rowptr[row_start+i];
            index_type local_k = this->partition[p].rowptr[i];
            for (int k = 0; k < nvals; ++k)
            {
               this->partition[p].colidx[local_k+k] = M.colidx[global_k+k];
               this->partition[p].values[local_k+k] = M.values[global_k+k];
               //++nops;
            }

            this->partition[p].rowptr[i+1] = this->partition[p].rowptr[i] + nvals;
         }

         if (__DEBUG) printf("partition[%d]: row_start=%d, row_stop=%d, num_rows=%d, num_columns=%d, num_values=%d\n", p, row_start, row_stop, this->partition[p].num_rows, this->partition[p].num_columns, this->partition[p].num_values);
      }
      //printf("ncopied = %d\n", nops);

#ifdef USE_LIBNUMA
      if (parallel::numa_params.is_initialized)
         numa_set_membind(&default_mask);
#endif
   }

   template <typename Format>
   void _copy(const MatrixType &M, Format format)
   {
      fprintf(stderr,"format %s not supported here %s %d\n", splib::format_name(Format()).c_str(), __FILE__, __LINE__);
      exit(-1);
   }

   //partitioned_matrix (const MatrixType &M, int _num_partitions = parallel::numa_params.num_nodes)
   partitioned_matrix (const MatrixType &M, int _num_partitions = DefaultNumberPartitions)
      :
      Parent (M.num_rows, M.num_columns, M.num_values),
      num_partitions(_num_partitions), partition(NULL), offset(_num_partitions+1)
   {
      //printf("... constructing partitioned matrix ... %s\n", typeid(allocator_type).name().c_str());
      if (__DEBUG) std::cout << "... constructing (copy) partitioned matrix ... " << typeid(*this).name() << " " << num_partitions << std::endl;

      //partition = (matrix_type *) malloc(sizeof(matrix_type)*num_partitions);
      partition = new matrix_type [num_partitions];
      if (partition == NULL)
      {
         fprintf(stderr,"partitioned_matrix allocation failure\n");
         exit(-1);
      }

      parallel::partitioner(&this->offset[0], this->num_partitions, 0, this->num_rows);

      this->_copy(M, typename matrix_type::format_type());
   }

   partitioned_matrix (void)
      :
      Parent (),
      num_partitions(0), partition(NULL), offset()
   {}

   ~partitioned_matrix()
   {
      this->offset.clear();

      if (partition)
         delete [] partition;
         //free(partition);
   }

   //void copy (const MatrixType &M, int _num_partitions = parallel::numa_params.num_nodes)
   void copy (const MatrixType &M, int _num_partitions = DefaultNumberPartitions)
   {
      this->~partitioned_matrix();

      this->num_rows    = M.num_rows;
      this->num_columns = M.num_columns;
      this->num_values  = M.num_values;

      this->num_partitions = _num_partitions;

      if (__DEBUG) std::cout << "... copying partitioned matrix ... " << typeid(*this).name() << " " << num_partitions << std::endl;

      this->partition = new matrix_type [num_partitions];
      if (this->partition == NULL)
      {
         fprintf(stderr,"partitioned_matrix allocation failure %s %d\n", __FILE__,__LINE__);
         exit(-1);
      }

      this->offset.resize(num_partitions+1);

      parallel::partitioner(&this->offset[0], this->num_partitions, 0, this->num_rows);

      this->_copy(M, typename matrix_type::format_type());
   }
};

struct thread_info_type
{
   int num_threads;
   int thread_id;
   int partition_id;
   int partition_rank;
   int threads_per_partition;

   thread_info_type (void) : num_threads(1), thread_id(0), partition_id(0), partition_rank(0), threads_per_partition(0) {}

   thread_info_type (int _num_threads, int _thread_id, int _partition_id, int _partition_rank) :
      num_threads(_num_threads),
      thread_id(_thread_id),
      partition_id(_partition_id),
      partition_rank(_partition_rank)
   {}

   ~thread_info_type()
   {
      if (__DEBUG) { printf("destroyed thread_info_type() %x\n", this); }
   }
};

//static std::vector< thread_info_type > thread_info_objects;
//static const thread_info_type serial_thread_object;

namespace details
{

   template <typename Matrix, typename Vector1, typename Vector2>
   void serial_spmv (const Matrix& A, const Vector1& x, Vector2& y, const int row_start, const int row_stop, splib::csr_format)
   {
      typedef typename Vector2::value_type value_type;
      #pragma ivdep
      for (int i = row_start; i < row_stop; ++i)
      {
         value_type sum(0);
         const int k_stop = A.rowptr[i+1];
         for (int k = A.rowptr[i]; k < k_stop; ++k)
         {
            sum += A.values[k] * x[ A.colidx[k] ];
         }

         y[i] = sum;
      }
   }
   template <typename Matrix, typename Vector1, typename Vector2, typename Format>
   void serial_spmv (const Matrix& A, const Vector1& x, Vector2& y, const int row_start, const int row_stop, Format format)
   {
      fprintf(stderr,"format %s not supported here %s %d\n", splib::format_name(Format()).c_str(), __FILE__, __LINE__);
      exit(-1);
   }
   template <typename Matrix, typename Vector1, typename Vector2, typename Vector3>
   void serial_gemv (const typename Vector3::value_type alpha,
              const Matrix&	A,
              const Vector1&	x,
              const typename Vector3::value_type beta,
              const Vector2&	y,
                    Vector3&	z,
              const int row_start,
              const int row_stop,
              splib::csr_format)
   {
      typedef typename Matrix::index_type index_type;
      typedef typename Vector3::value_type value_type;

      const int num_rows    = A.num_rows;
      const int num_columns = A.num_columns;
      const int num_values  = A.num_values;

      #pragma ivdep
      for (int i = row_start; i < row_stop; ++i)
      {
         value_type sum(0);

         // A*x for this row ...
         const int k_stop  = A.rowptr[i+1];
         for (int k = A.rowptr[i]; k < k_stop; ++k)
            sum += (A.values[k] * x[ A.colidx[k] ]);

         // z := alpha*A*x + beta*y
         if (alpha == value_type(1) && beta == value_type(1))
            z[i] = sum + y[i];
         else if (alpha == value_type(1) && beta == value_type(-1))
            z[i] = sum - y[i];
         else if (alpha == value_type(-1) && beta == value_type(1))
            z[i] = y[i] - sum;
         else if (alpha == value_type(1))
            z[i] = sum + beta*y[i];
         else if (beta == value_type(1))
            z[i] = alpha*sum + y[i];
         else
            z[i] = alpha*sum + beta*y[i];
      }

      return;
   }
   template <typename Matrix, typename Vector1, typename Vector2, typename Vector3, typename Format>
   void serial_gemv (const typename Vector3::value_type alpha,
              const Matrix&	A,
              const Vector1&	x,
              const typename Vector3::value_type beta,
              const Vector2&	y,
                    Vector3&	z,
              const int row_start,
              const int row_stop,
              const Format format)
   {
      fprintf(stderr,"serial_gemv: format %s not supported here %s %d\n", splib::format_name(Format()).c_str(), __FILE__, __LINE__);
      exit(-1);
   }

   // SpMV ... x has global index space!!!
   template <typename Matrix, typename Vector1, typename Vector2, typename ThreadInfo>
   void multiply (const Matrix& A,
                  const Vector1& x,
                        Vector2& y,
                  const ThreadInfo& thread_info,
                          splib::partitioned_matrix_format,
                          splib::vector_format,
                          splib::partitioned_vector_format)
   {
      typedef typename Vector2::value_type value_type;

      assert(A.num_partitions == y.num_partitions);

      //const int thread_id = parallel::get_thread_num();
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         // make sure the global x is finished writing.
         #pragma omp barrier

         // Nested parallelism ... let the kernel handle the parallel_for()
         if (num_threads == A.num_partitions)
         {
            splib::multiply( A.partition[thread_id], x, y.partition[thread_id]);
            //return;
         }
         // Persistent parallelism
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            const int row_start = A.offset[p];
            const int row_stop  = A.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Matrix::matrix_type& Ap = A.partition[p];
            const typename Vector2::vector_type& yp = y.partition[p];

            details::serial_spmv (Ap, x, yp, i_start, i_stop, typename Matrix::matrix_type::format_type());
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         for (int p = 0; p < A.num_partitions; ++p)
         {
            assert(A.offset[p] == y.offset[p]);

            splib::multiply (A.partition[p], x, y.partition[p]);
         }
   }

   template <typename Vector1, typename Vector2, typename ValueType, typename ThreadInfo>
   ValueType dot (const Vector1& x,
                  const Vector2& y,
                        ValueType& sum,
                  const ThreadInfo& thread_info,
                              splib::partitioned_vector_format,
                              splib::partitioned_vector_format)
   {
      //typedef typename Vector::value_type value_type;
      typedef ValueType value_type;

      assert(x.num_partitions == y.num_partitions);

      #pragma omp single
      {
         sum = ValueType(0);
      }
      // implicit barrier

      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         ValueType local_sum(0);

         if (num_threads == x.num_partitions) // nested ..
         {
            local_sum = splib::spblas::dot(x.partition[thread_id], y.partition[thread_id]);
         }
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            assert(x.offset[p] == y.offset[p]);

            const int row_start = x.offset[p];
            const int row_stop  = x.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Vector1::vector_type& xp = x.partition[p];
            const typename Vector2::vector_type& yp = y.partition[p];

            #pragma ivdep
            for (int i = i_start; i < i_stop; ++i)
               local_sum += (xp[i] * yp[i]);
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }

         #pragma omp critical
         {
            sum += local_sum;
         }
         #pragma omp barrier
      }
      else
         for (int p = 0; p < x.num_partitions; ++p)
         {
            assert(x.offset[p] == y.offset[p]);

            sum += splib::spblas::dot( x.partition[p], y.partition[p]);
         }

      return sum;
   }
   template <typename Vector1, typename Vector2, typename ValueType, typename ThreadInfo>
   ValueType dot (const Vector1& x,
                  const Vector2& y,
                        ValueType& sum,
                  const ThreadInfo& thread_info,
                             splib::vector_format,
                             splib::partitioned_vector_format)
   {
      #pragma omp single
      { sum = 0; }
      // implicit barrier

      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         ValueType local_sum(0);

         const int thread_id = parallel::get_thread_num();

         if (num_threads == y.num_partitions) // nested ..
         {
            int row_start = y.offset[thread_id];
            int row_stop  = y.offset[thread_id+1];
            int num_rows  = row_stop - row_start;

            Vector1 xp( &x[row_start], num_rows);

            local_sum = splib::spblas::dot(xp, y.partition[thread_id]);
         }
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            const int row_start = y.offset[p];
            const int row_stop  = y.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Vector2::vector_type& yp = y.partition[p];

            #pragma ivdep
            for (int i = i_start; i < i_stop; ++i)
               local_sum += x[row_start + i] * yp[i];
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }

         #pragma omp critical
         { sum += local_sum; }
         #pragma omp barrier
      }
      else
         for (int p = 0; p < y.num_partitions; ++p)
         {
            int row_start = y.offset[p];
            int row_stop  = y.offset[p+1];
            int num_rows  = row_stop - row_start;

            Vector1 xp( &x[row_start], num_rows);

            sum += splib::spblas::dot( xp, y.partition[p]);
         }

      return sum;
   }
   template <typename Vector1, typename Vector2, typename ValueType, typename ThreadInfo>
   ValueType dot (const Vector1& x,
                  const Vector2& y,
                        ValueType& sum,
                  const ThreadInfo& thread_info,
                             splib::partitioned_vector_format,
                             splib::vector_format)
   {
      return dot(y,x,sum,thread_info);
   }

   template <typename Value1, typename Vector1, typename Vector2, typename Vector3>
   void axpby (const Value1 alpha,
               const Value1 beta,
               const Vector1& x,
               const Vector2& y,
                     Vector3& z,
               const thread_info_type& thread_info,
                    splib::partitioned_vector_format, // x
                    splib::partitioned_vector_format, // y
                    splib::partitioned_vector_format) // z
   {
      assert(x.num_partitions == y.num_partitions);
      assert(z.num_partitions == y.num_partitions);

      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         if (num_threads == x.num_partitions) // nested ..
         {
            splib::spblas::axpby (alpha, beta, x.partition[thread_id], y.partition[thread_id], z.partition[thread_id]);
         }
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            assert(x.offset[p] == y.offset[p]);

            const int row_start = x.offset[p];
            const int row_stop  = x.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Vector1::vector_type& xp = x.partition[p];
            const typename Vector2::vector_type& yp = y.partition[p];
                  typename Vector3::vector_type& zp = z.partition[p];

            //#pragma omp for
            //for (int i = 0; i < num_rows; ++i)
            #pragma ivdep
            for (int i = i_start; i < i_stop; ++i)
            {
               if (alpha == Value1(1) && beta == Value1(1))
                  zp[i] = xp[i] + yp[i];
               else if (alpha == Value1(1) && beta == Value1(-1))
                  zp[i] = xp[i] - yp[i];
               else if (alpha == Value1(-1) && beta == Value1(1))
                  zp[i] = yp[i] - xp[i];
               else if (alpha == Value1(1))
                  zp[i] = xp[i] + beta*yp[i];
               else if (beta == Value1(1))
                  zp[i] = alpha*xp[i] + yp[i];
               else
                  zp[i] = alpha*xp[i] + beta*yp[i];
            }
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         for (int p = 0; p < x.num_partitions; ++p)
         {
            assert(x.offset[p] == y.offset[p]);
            assert(z.offset[p] == y.offset[p]);

            splib::spblas::axpby(alpha, beta, x.partition[p], y.partition[p], z.partition[p]);
         }
   }
   template <typename Value1, typename Vector1, typename Vector2, typename Vector3>
   void axpby (const Value1 alpha,
               const Value1 beta,
               const Vector1& x,
               const Vector2& y,
                     Vector3& z,
               const thread_info_type& thread_info,
                    splib::partitioned_vector_format, // x
                    splib::vector_format, // y
                    splib::vector_format) // z
   {
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         if (num_threads == x.num_partitions) // nested ..
         {
            int row_start = x.offset[thread_id];
            int row_stop  = x.offset[thread_id+1];
            int num_rows  = row_stop - row_start;

            Vector2 yp(&y[row_start], num_rows);
            Vector3 zp(&z[row_start], num_rows);

            splib::spblas::axpby (alpha, beta, x.partition[thread_id], yp, zp);
         }
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            const int row_start = x.offset[p];
            const int row_stop  = x.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Vector1::vector_type& xp = x.partition[p];

            //#pragma omp for
            //for (int i = 0; i < num_rows; ++i)
            #pragma ivdep
            for (int i = i_start; i < i_stop; ++i)
            {
               if (alpha == Value1(1) && beta == Value1(1))
                  z[row_start+i] = xp[i] + y[row_start+i];
               else if (alpha == Value1(1) && beta == Value1(-1))
                  z[row_start+i] = xp[i] - y[row_start+i];
               else if (alpha == Value1(-1) && beta == Value1(1))
                  z[row_start+i] = y[row_start+i] - xp[i];
               else if (alpha == Value1(1))
                  z[row_start+i] = xp[i] + beta*y[row_start+i];
               else if (beta == Value1(1))
                  z[row_start+i] = alpha*xp[i] + y[row_start+i];
               else
                  z[row_start+i] = alpha*xp[i] + beta*y[row_start+i];
            }
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         for (int p = 0; p < x.num_partitions; ++p)
         {
            int row_start = x.offset[p];
            int row_stop  = x.offset[p+1];
            int num_rows  = row_stop - row_start;

            Vector2 yp(&y[0] + row_start, num_rows);
            Vector3 zp(&z[0] + row_start, num_rows);

            splib::spblas::axpby (alpha, beta, x.partition[p], yp, zp);
         }
   }
   template <typename Value1, typename Vector1, typename Vector2, typename Vector3>
   void axpby (const Value1 alpha,
               const Value1 beta,
               const Vector1& x,
               const Vector2& y,
                     Vector3& z,
               const thread_info_type& thread_info,
                    splib::vector_format, // x
                    splib::vector_format, // y
                    splib::vector_format) // z
   {
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         const int num_rows = x.size();

         if (num_threads == parallel::num_threads) // persistent threads ...
         {
            #pragma ivdep
            #pragma omp for
            for (int i = 0; i < num_rows; ++i)
            {
               if (alpha == Value1(1) && beta == Value1(1))
                  z[i] = x[i] + y[i];
               else if (alpha == Value1(1) && beta == Value1(-1))
                  z[i] = x[i] - y[i];
               else if (alpha == Value1(-1) && beta == Value1(1))
                  z[i] = y[i] - x[i];
               else if (alpha == Value1(1))
                  z[i] = x[i] + beta*y[i];
               else if (beta == Value1(1))
                  z[i] = alpha*x[i] + y[i];
               else
                  z[i] = alpha*x[i] + beta*y[i];
            }
         }
         else // nested ..
         {
            int row_start, row_stop;
            parallel::partitioner (row_start, row_stop, thread_id, num_threads, 0, num_rows);

            int _num_rows  = row_stop - row_start;

            Vector1 xp(&x[row_start], _num_rows);
            Vector2 yp(&y[row_start], _num_rows);
            Vector3 zp(&z[row_start], _num_rows);

            splib::spblas::axpby (alpha, beta, xp, yp, zp);
         }
         //else
         //   { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         splib::spblas::axpby (alpha, beta, x, y, z);
   }

   template <typename Matrix, typename Vector1, typename Vector2, typename Vector3>
   void gemv (const typename Vector3::value_type alpha,
              const Matrix& A,
              const Vector1& x,
              const typename Vector3::value_type beta,
              const Vector2& y,
                    Vector3& z,
              const thread_info_type& thread_info,
                    splib::partitioned_matrix_format, // A
                    splib::vector_format, // x
                    splib::partitioned_vector_format, // y
                    splib::partitioned_vector_format) // z
   {
      //details::multiply (A, x, z, thread_info);
      //details::axpby (alpha, beta, z, y, z, thread_info);
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         // make sure the global x is finished writing.
         //#pragma omp barrier

         // Nested parallelism ... let the kernel handle the parallel_for()
         if (num_threads == A.num_partitions)
         {
            splib::gemv( alpha, A.partition[thread_id], x, beta, y.partition[thread_id], z.partition[thread_id]);
            //return;
         }
         // Persistent parallelism
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            const int row_start = A.offset[p];
            const int row_stop  = A.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Matrix::matrix_type& Ap = A.partition[p];
            const typename Vector2::vector_type& yp = y.partition[p];
            const typename Vector3::vector_type& zp = z.partition[p];

            details::serial_gemv (alpha, Ap, x, beta, yp, zp, i_start, i_stop, typename Matrix::matrix_type::format_type());
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         for (int p = 0; p < A.num_partitions; ++p)
         {
            assert(A.offset[p] == y.offset[p]);

            //splib::multiply (A.partition[p], x, y.partition[p]);
            splib::gemv( alpha, A.partition[p], x, beta, y.partition[p], z.partition[p]);
         }
   }
   template <typename Matrix, typename Vector1, typename Vector2, typename Vector3>
   void gemv (const typename Vector3::value_type alpha,
              const Matrix& A,
              const Vector1& x,
              const typename Vector3::value_type beta,
              const Vector2& y,
                    Vector3& z,
              const thread_info_type& thread_info,
                    splib::partitioned_matrix_format, // A
                    splib::vector_format, // x
                    splib::partitioned_vector_format, // y
                    splib::vector_format) // z
   {
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         // make sure the global x is finished writing.
         //#pragma omp barrier

         // Nested parallelism ... let the kernel handle the parallel_for()
         if (num_threads == A.num_partitions)
         {
            const int row_start = A.offset[thread_id];
            const int row_stop  = A.offset[thread_id+1];
            const int num_rows  = row_stop - row_start;

            Vector3 zp( &z[row_start], num_rows);

            splib::gemv( alpha, A.partition[thread_id], x, beta, y.partition[thread_id], zp);
            //return;
         }
         // Persistent parallelism
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            const int row_start = A.offset[p];
            const int row_stop  = A.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Matrix::matrix_type& Ap = A.partition[p];
            const typename Vector2::vector_type& yp = y.partition[p];
            Vector3 zp( &z[row_start], num_rows);

            details::serial_gemv (alpha, Ap, x, beta, yp, zp, i_start, i_stop, typename Matrix::matrix_type::format_type());
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         for (int p = 0; p < A.num_partitions; ++p)
         {
            assert(A.offset[p] == y.offset[p]);

            const int row_start = A.offset[thread_id];
            const int row_stop  = A.offset[thread_id+1];
            const int num_rows  = row_stop - row_start;

            Vector3 zp( &z[row_start], num_rows);

            splib::gemv( alpha, A.partition[p], x, beta, y.partition[p], zp);
         }
   }
   template <typename Matrix, typename Vector1, typename Vector2, typename Vector3>
   void gemv (const typename Vector3::value_type alpha,
              const Matrix& A,
              const Vector1& x,
              const typename Vector3::value_type beta,
              const Vector2& y,
                    Vector3& z,
              const thread_info_type& thread_info,
                    splib::partitioned_matrix_format, // A
                    splib::vector_format, // x
                    splib::vector_format, // y
                    splib::vector_format  // z
                  )
   {
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         // make sure the global x is finished writing.
         //#pragma omp barrier

         // Nested parallelism ... let the kernel handle the parallel_for()
         if (num_threads == A.num_partitions)
         {
            const int row_start = A.offset[thread_id];
            const int row_stop  = A.offset[thread_id+1];
            const int num_rows  = row_stop - row_start;

            Vector2 yp(&y[row_start], num_rows);
            Vector3 zp(&z[row_start], num_rows);

            splib::gemv( alpha, A.partition[thread_id], x, beta, yp, zp);
         }
         // Persistent parallelism
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            const int row_start = A.offset[p];
            const int row_stop  = A.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Matrix::matrix_type& Ap = A.partition[p];

            Vector2 yp(&y[row_start], num_rows);
            Vector3 zp(&z[row_start], num_rows);

            details::serial_gemv (alpha, Ap, x, beta, yp, zp, i_start, i_stop, typename Matrix::matrix_type::format_type());
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         for (int p = 0; p < A.num_partitions; ++p)
         {
            //assert(A.offset[p] == y.offset[p]);

            const int row_start = A.offset[p];
            const int row_stop  = A.offset[p+1];
            const int num_rows  = row_stop - row_start;

            Vector2 yp(&y[row_start], num_rows);
            Vector3 zp(&z[row_start], num_rows);

            splib::gemv( alpha, A.partition[p], x, beta, yp, zp);
         }
   }

   template <typename Vector1, typename Vector2, typename Vector3, typename ThreadInfo>
   void xmy (const Vector1& x,
             const Vector2& y,
                   Vector3& z,
             const ThreadInfo& thread_info,
                    splib::partitioned_vector_format, // x
                    splib::vector_format, // y
                    splib::vector_format) // z
   {
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         if (num_threads == x.num_partitions) // nested ..
         {
            int row_start = x.offset[thread_id];
            int row_stop  = x.offset[thread_id+1];
            int num_rows  = row_stop - row_start;

            Vector2 yp(&y[row_start], num_rows);
            Vector3 zp(&z[row_start], num_rows);

            splib::spblas::xmy (x.partition[thread_id], yp, zp);
         }
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            const int row_start = x.offset[p];
            const int row_stop  = x.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Vector1::vector_type& xp = x.partition[p];

            #pragma ivdep
            for (int i = i_start; i < i_stop; ++i)
            {
               z[row_start+i] = xp[i] * y[row_start+i];
            }
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         for (int p = 0; p < x.num_partitions; ++p)
         {
            int row_start = x.offset[p];
            int row_stop  = x.offset[p+1];
            int num_rows  = row_stop - row_start;

            Vector2 yp(&y[0] + row_start, num_rows);
            Vector3 zp(&z[0] + row_start, num_rows);

            splib::spblas::xmy (x.partition[p], yp, zp);
         }
   }
   template <typename Vector1, typename Vector2, typename Vector3, typename ThreadInfo>
   void xmy (const Vector1& x,
             const Vector2& y,
                   Vector3& z,
             const ThreadInfo& thread_info,
                    splib::vector_format, // x
                    splib::vector_format, // y
                    splib::vector_format) // z
   {
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      const int num_rows = x.size();

      if (num_threads > 1)
      {
         if (num_threads == parallel::num_threads)
         {
            #pragma omp for
            for (int i = 0; i < num_rows; ++i)
            {
               z[i] = x[i] * y[i];
            }
         }
         else
            if (EnableNestedParallelism)// and num_threads == thread_info.num_partitions) // nested ..
            {
               int i_start, i_stop;
               parallel::partitioner (i_start, i_stop, thread_id, num_threads, 0, num_rows);

               const int nrows = i_stop - i_start;

               Vector1 xp(&x[i_start], nrows);
               Vector2 yp(&y[i_start], nrows);
               Vector3 zp(&z[i_start], nrows);

               #pragma omp for
               for (int i = 0; i < nrows; ++i)
               {
                  zp[i] = xp[i] * yp[i];
               }
            }
            else
               { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         splib::spblas::xmy (x, y, z);
   }

   template <typename Vector1, typename Vector2, typename ThreadInfo>
   void copy (const Vector1& x,
                    Vector2& y,
              const ThreadInfo &thread_info,
                    splib::partitioned_vector_format, // x
                    splib::vector_format) // y
   {
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      if (num_threads > 1)
      {
         if (num_threads == x.num_partitions) // nested ..
         {
            int row_start = x.offset[thread_id];
            int row_stop  = x.offset[thread_id+1];
            int num_rows  = row_stop - row_start;

            Vector2 yp(&y[row_start], num_rows);

            splib::spblas::copy (x.partition[thread_id], yp);
         }
         else if (num_threads == parallel::num_threads)
         {
            const int& p = thread_info.partition_id;
            const int& tid = thread_info.partition_rank;
            const int& threads_per_partition = thread_info.threads_per_partition;

            const int row_start = x.offset[p];
            const int row_stop  = x.offset[p+1];
            const int num_rows  = row_stop - row_start;

            int i_start, i_stop;
            parallel::partitioner (i_start, i_stop, tid, threads_per_partition, 0, num_rows);

            const typename Vector1::vector_type& xp = x.partition[p];

            #pragma ivdep
            for (int i = i_start; i < i_stop; ++i)
            {
               y[row_start+i] = xp[i];
            }
         }
         else
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         for (int p = 0; p < x.num_partitions; ++p)
         {
            int row_start = x.offset[p];
            int row_stop  = x.offset[p+1];
            int num_rows  = row_stop - row_start;

            Vector2 yp(&y[0] + row_start, num_rows);

            splib::spblas::copy (x.partition[p], yp);
         }
   }
   template <typename Vector1, typename Vector2, typename ThreadInfo>
   void copy (const Vector1& x,
                    Vector2& y,
              const ThreadInfo &thread_info,
                    splib::vector_format, // x
                    splib::vector_format) // y
   {
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      const int len = x.size();

      if (num_threads > 1)
      {
#ifdef _OPENMP
         if (not(omp_in_parallel()))
            { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
#endif
         #pragma omp for
         for (int i = 0; i < len; ++i)
         {
            y[i] = x[i];
         }
      }
      else
         splib::spblas::copy (x, y);
   }

   template <typename ValueType, typename Vector1, typename Vector2, typename ThreadInfo>
   void scale (const ValueType& alpha,
               const Vector1& x,
                     Vector2& y,
               const ThreadInfo &thread_info,
                     splib::vector_format, // x
                     splib::vector_format) // y
   {
      const int& thread_id = thread_info.thread_id;
      const int& num_threads = thread_info.num_threads;

      const int num_rows = x.size();

      if (num_threads > 1)
      {
         if (num_threads == parallel::num_threads)
         {
            #pragma omp for
            for (int i = 0; i < num_rows; ++i)
            {
               y[i] = alpha * x[i];
            }
         }
         else
            if (EnableNestedParallelism) // num_threads should == the # of nodes
            {
               int i_start, i_stop;
               parallel::partitioner (i_start, i_stop, thread_id, num_threads, 0, num_rows);

               const int nrows = i_stop - i_start;

               Vector1 xp(&x[i_start], nrows);
               Vector2 yp(&y[i_start], nrows);

               #pragma omp for
               for (int i = 0; i < nrows; ++i)
               {
                  yp[i] = alpha * xp[i];
               }
            }
            else
               { fprintf(stderr,"this should not happen %s %d\n", __FILE__, __LINE__); exit(-1); }
      }
      else
         splib::spblas::scale (alpha, x, y);
   }

} // namespace details

   template <typename Matrix, typename Vector1, typename Vector2, typename Vector3>
   void gemv (const typename Vector3::value_type alpha,
              const Matrix& A,
              const Vector1& x,
              const typename Vector3::value_type beta,
              const Vector2& y,
                    Vector3& z,
              const thread_info_type& thread_info)
   {
      details::gemv(alpha, A, x, beta, y, z, thread_info,
                       typename Matrix::format_type(),
                       typename Vector1::format_type(),
                       typename Vector2::format_type(),
                       typename Vector3::format_type());
   }

   template <typename Matrix, typename Vector1, typename Vector2, typename ThreadInfo>
   void multiply (const Matrix& A,
                  const Vector1& x,
                        Vector2& y,
                  const ThreadInfo& thread_info)
   {
      details::multiply(A, x, y, thread_info,
                       typename Matrix::format_type(),
                       typename Vector1::format_type(),
                       typename Vector2::format_type());
   }

   template <typename Vector1, typename Vector2, typename ValueType, typename ThreadInfo>
   ValueType dot (const Vector1& x, const Vector2& y, ValueType& sum, const ThreadInfo& thread_info)
   {
      return details::dot(x, y, sum, thread_info,
                             typename Vector1::format_type(),
                             typename Vector2::format_type());
   }

   template <typename Vector, typename ValueType, typename ThreadInfo>
   ValueType norm2 (const Vector& v, ValueType& sum, const ThreadInfo& thread_info)
   {
      dot(v, v, sum, thread_info);

      #pragma omp single
      {
         //printf("dot sum = %f %f\n", sum, std::sqrt(sum));
         sum = std::sqrt(sum);
      }
      // implicit barrier

      return sum;
   }

   template <typename Vector1, typename Vector2, typename Vector3, typename ThreadInfo>
   void xmy (const Vector1& x, const Vector2& y, Vector3& z, const ThreadInfo& thread_info)
   {
      details::xmy(x, y, z, thread_info,
                             typename Vector1::format_type(),
                             typename Vector2::format_type(),
                             typename Vector3::format_type());
   }

   template <typename ValueType, typename Vector1, typename Vector2, typename Vector3>
   void axpby (const ValueType alpha,
               const ValueType beta,
               const Vector1& x,
               const Vector2& y,
                     Vector3& z,
               const thread_info_type& thread_info)
   {
      details::axpby (alpha, beta, x, y, z, thread_info,
                             typename Vector1::format_type(),
                             typename Vector2::format_type(),
                             typename Vector3::format_type());
   }

   template <typename Vector1, typename Vector2, typename ThreadInfo>
   void copy (const Vector1& x, Vector2& y, const ThreadInfo& thread_info)
   {
      details::copy (x, y, thread_info,
                             typename Vector1::format_type(),
                             typename Vector2::format_type());
   }
   template <typename Vector1, typename Vector2>
   void copy (const partitioned_vector<Vector1>& x, Vector2& y)
   {
      bool in_parallel(false);

#ifdef _OPENMP
      in_parallel = omp_in_parallel();
#endif
      #pragma omp single
      for (int p = 0; p < x.num_partitions; ++p)
      {
         int row_start = x.offset[p];
         int row_stop  = x.offset[p+1];
         int num_rows  = row_stop - row_start;

         Vector2 yp(&y[0] + row_start, num_rows);

         splib::spblas::copy (x.partition[p], yp);
      }
   }

   template <typename Vector1, typename Vector2, typename ThreadInfo>
   void scale (const typename Vector1::value_type& alpha, const Vector1& x, Vector2& y, const ThreadInfo& thread_info)
   {
      details::scale (alpha, x, y, thread_info,
                             typename Vector1::format_type(),
                             typename Vector2::format_type());
   }

/*   template <typename Vector1, typename Vector2>
   void copy (const Vector1& x, Vector2& y)
   {
      parallel::thread_info_type thread_info;

#ifdef _OPENMP
      thread_info.num_threads = omp_get_num_threads();
      thread_info.thread_id = omp_get_thread_num();
#endif
      if (EnableNestedParallelism)
         thread_info.partition_id = thread_info.thread_id;
      else
      {
         thread_info.threads_per_partition = thread_info.num_threads / x.num_partitions;
         thread_info.partition_id = thread_info.thread_id / thread_info.threads_per_partition;
         thread_info.partition_rank = thread_info.thread_id - thread_info.partition_id * thread_info.threads_per_partition;
      }
      printf("copy: thread_info %d %d %d %d\n", thread_info.thread_id, thread_info.partition_id, thread_info.partition_rank, thread_info.num_threads);

      copy(x,y,thread_info);
   }*/


} // namespace parallel

#endif
