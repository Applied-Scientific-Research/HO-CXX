#include <cstdio>
#include <cstdlib>
#include <type_traits>

#include "aplles_interface.h"

template <typename T>
T* allocate ( const size_t nelems )
{
   T *ptr = (T*) malloc( nelems * sizeof(T) );
   if (ptr == NULL)
   {
      fprintf(stderr,"Allocation error for %lu elements\n", nelems);
      return NULL;
   }
   else
      return ptr;
}

template <typename T>
int allocate ( const size_t nelems, T* &ptr )
{
   ptr = allocate<T>( nelems );
   if (ptr == NULL)
   {
      fprintf(stderr,"Allocation error for %lu elements at %d %s\n", nelems, __LINE__, __FILE__);
      exit(-1);
      return 1;
   }
   else
      return 0;
}

template <typename T>
void deallocate ( T* &ptr )
{
   if (ptr)
   {
      free(ptr);
      ptr = NULL;
   }
}

template <typename T>
int reallocate ( const size_t nelems, T* &ptr )
{
   T *new_ptr = (T*) realloc( (void *) ptr, nelems * sizeof(T) );
   if ( new_ptr == NULL )
   {
      fprintf(stderr,"Reallocation error for %lu elements at %d %s\n", nelems, __LINE__, __FILE__);
      deallocate( ptr );
      exit(-1);
      return 1;
   }
   else
   {
      ptr = new_ptr;
      return 0;
   }
}

template <typename T>
struct DynamicArrayType
{
   typedef T ValueType;

   ValueType *m_ptr;
   bool       m_owns_data; /// Is this a reference type that points to another allocation?
   size_t     m_nelems;  /// How many elements?

   explicit DynamicArrayType ( ValueType *ptr, const size_t nelems ) : m_ptr( ptr ), m_nelems(nelems), m_owns_data(false)
   {
      fprintf(stderr,"Instantiated DynamicArrayType %p with existing array %p\n", this, this->m_ptr);
   }

   DynamicArrayType (void) : m_ptr(NULL), m_nelems(0), m_owns_data(true)
   {
      fprintf(stderr,"Instantiated empty DynamicArrayType %p\n", this);
   }

   DynamicArrayType ( const size_t nelems) : m_ptr(NULL), m_nelems(0), m_owns_data(true)
   {
      allocate( nelems, this->m_ptr );
      fprintf(stderr,"Instantiated DynamicArrayType %p with %lu elements %p\n", this, this->m_nelems, this->m_ptr);
   }

   ~DynamicArrayType (void)
   {
      this->clearData();
   }

   void clearData (void )
   {
      if ( this->m_owns_data and this->m_nelems > 0 and this->m_ptr )
      {
         fprintf(stderr,"Clearing object's internal data %p %p %lu\n", this, this->m_ptr, this->m_nelems);
         deallocate( this->m_ptr );
         this->m_nelems = 0;
      }
   }

   void setData ( ValueType *ptr, const size_t nelems )
   {
      this->clearData();
      this->m_nelems = nelems;
      this->m_owns_data = false;
      this->m_ptr = ptr;
      fprintf(stderr,"DynamicArrayType::setData %p with existing array %p %lu\n", this, this->m_ptr, this->m_nelems);
   }

   void resize ( const size_t nelems )
   {
      if ( this->m_owns_data == false )
      {
          fprintf(stderr,"Error: attempting to resize a referenced array %p\n", this);
          exit(-1);
      }
      else
      {
          reallocate( nelems, this->m_ptr );
          this->m_nelems = nelems;
      }
   }

         T& operator() (const size_t i)       { return this->m_ptr[i]; }
   const T& operator() (const size_t i) const { return this->m_ptr[i]; }
};

template < typename T, int ... dims >
struct StaticArrayProperties;

template < typename T, int _I0 >
struct StaticArrayProperties< T[_I0] >
{
   enum { rank = 1 };
   typedef T ValueType;
   enum { I0 = _I0 };
   enum { I1 = 0 };
};

template < typename T, int _I0, int _I1 >
struct StaticArrayProperties< T[_I1][_I0] >
{
   enum { rank = 2 };
   typedef T ValueType;
   enum { I0 = _I0 };
   enum { I1 = _I1 };
};

template < class DataType >
struct StaticArrayType
{
   typedef StaticArrayProperties<DataType> properties;

   enum { rank = properties::rank };
   enum { I0   = properties::I0 };
   enum { I1   = properties::I1 };

   typedef typename properties::ValueType ValueType;

   DataType m_data;

   /// rank 2
   template <typename Int>
      typename std::enable_if< std::is_integral<Int>::value and rank == 2, ValueType& >::type
   operator()( const Int i, const Int j )
   {
      return this->m_data[j][i];
   }

   template <typename Int>
     const typename std::enable_if< std::is_integral<Int>::value and rank == 2, ValueType& >::type
   operator()( const Int i, const Int j ) const
   {
      return this->m_data[j][i];
   }

   /// rank 1
   template <typename Int>
      typename std::enable_if< std::is_integral<Int>::value and rank == 1, ValueType& >::type
   operator()( const Int i )
   {
      return this->m_data[i];
   }

   template <typename Int>
     const typename std::enable_if< std::is_integral<Int>::value and rank == 1, ValueType& >::type
   operator()( const Int i ) const
   {
      return this->m_data[i];
   }

   /// rank 1 ... to match normal access format.
   template <typename Int>
      typename std::enable_if< std::is_integral<Int>::value and rank == 1, ValueType& >::type
   operator[]( const Int i )
   {
      return this->m_data[i];
   }

   template <typename Int>
     const typename std::enable_if< std::is_integral<Int>::value and rank == 1, ValueType& >::type
   operator[]( const Int i ) const
   {
      return this->m_data[i];
   }
 
};

template <int _Knod>
struct LaplacianDataType
{
   enum { Knod = _Knod };

   typedef StaticArrayType< double[Knod] > ElementNodeType;

   int Nel; /// Number of elements in mesh.
   DynamicArrayType< ElementNodeType > Vol_Jac; /// Volumetric Jacobian of the element.

   ElementNodeType wgt;

   int NelB; /// Number of boundary
   DynamicArrayType< int > BoundaryPointElementIDs;

   typedef StaticArrayType< double[Knod*Knod][Knod] > BoundarySourceType;
   DynamicArrayType< BoundarySourceType > BoundarySource;

   typedef StaticArrayType< double[Knod] > BoundaryValueType;
   DynamicArrayType< BoundaryValueType > BoundaryValues;

   /// APLLES handles
   APLLES_MatrixHandle_t MatrixHandle;
   APLLES_SolverHandle_t SolverHandle;
   APLLES_PreconHandle_t PreconHandle;
};

typedef LaplacianDataType<3> LaplacianType;
LaplacianType LaplacianData;

extern "C"
{

void setLaplacian( const int Knod, const int Nel,
                     double wgt[], double Vol_Jac[],
                     const int NelB, int BoundaryPointElementIDs1[], double BoundarySource[], double BoundaryValues[],
                     void *A_handle, void *S_handle, void *P_handle )
{
   printf("Knod: %d\n", Knod);
   printf("Nel: %d\n", Nel);
   printf("NelB: %d\n", NelB);

   if (Knod != 3) {
      fprintf(stderr,"Knod != 3 %d\n", Knod);
      exit(-1);
   }

   LaplacianData.Nel = Nel;

   for (int i(0); i < Knod; i++) {
      LaplacianData.wgt(i) = wgt[i];
      printf("wgt[%d]= %f\n", i, LaplacianData.wgt[i]);
   }

   LaplacianData.Vol_Jac.setData( (LaplacianType::ElementNodeType *) Vol_Jac, Nel );

   LaplacianData.NelB = NelB;

   LaplacianData.BoundarySource.setData( (LaplacianType::BoundarySourceType *) BoundarySource, NelB );
   LaplacianData.BoundaryValues.setData( (LaplacianType::BoundaryValueType  *) BoundaryValues, NelB );

   LaplacianData.BoundaryPointElementIDs.resize( NelB );
   for (int i(0); i < NelB; ++i)
      LaplacianData.BoundaryPointElementIDs(i) = BoundaryPointElementIDs1[i] - 1;

   LaplacianData.MatrixHandle = (APLLES_MatrixHandle_t) A_handle;
   LaplacianData.SolverHandle = (APLLES_SolverHandle_t) S_handle;
   LaplacianData.PreconHandle = (APLLES_PreconHandle_t) P_handle;

   return;
}

void getLaplacian( const double Vorticity[], double psi[] )
{
}

}
