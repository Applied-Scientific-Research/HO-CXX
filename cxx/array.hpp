#ifndef __array_hpp
#define __array_hpp

namespace HighOrderFEM
{

template <typename T>
struct getBaseValueType;

template <typename T>
struct DynamicArrayType
{
   typedef T ValueType;

   ValueType *m_ptr;
   bool       m_owns_data; /// Is this a reference type that points to another allocation?
   size_t     m_nelems;  /// How many elements?

   explicit DynamicArrayType ( ValueType *ptr, const size_t nelems ) : m_ptr( ptr ), m_nelems(nelems), m_owns_data(false)
   {
      //fprintf(stderr,"Instantiated DynamicArrayType %p with existing array %p\n", this, this->m_ptr);
   }

   DynamicArrayType (void) : m_ptr(NULL), m_nelems(0), m_owns_data(true)
   {
      //fprintf(stderr,"Instantiated empty DynamicArrayType %p\n", this);
   }

   DynamicArrayType ( const size_t nelems ) : m_ptr(NULL), m_nelems(0), m_owns_data(true)
   {
      allocate( nelems, this->m_ptr );
      //fprintf(stderr,"Instantiated DynamicArrayType %p with %lu elements %p\n", this, this->m_nelems, this->m_ptr);
   }

   ~DynamicArrayType (void)
   {
      this->clearData();
   }

   void clearData (void )
   {
      if ( this->m_owns_data and this->m_nelems > 0 and this->m_ptr )
      {
         //fprintf(stderr,"Clearing object's internal data %p %p %lu\n", this, this->m_ptr, this->m_nelems);
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
      //fprintf(stderr,"DynamicArrayType::setData %p with existing array %p %lu\n", this, this->m_ptr, this->m_nelems);
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

   ValueType *getPointer(void) { return this->m_ptr; }

   typename getBaseValueType<ValueType>::type* getRawPointer(void)
   {
      return (typename getBaseValueType<ValueType>::type* )(this->m_ptr);
   }

         T& operator() (const size_t i)       { return this->m_ptr[i]; }
   const T& operator() (const size_t i) const { return this->m_ptr[i]; }
         T& operator[] (const size_t i)       { return this->m_ptr[i]; }
   const T& operator[] (const size_t i) const { return this->m_ptr[i]; }
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

template <class ArrayType, int _Dim>
struct SliceType;

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
      typename std::enable_if< std::is_integral<Int>::value and rank == 2, ValueType >::type &
   operator()( const Int i, const Int j )
   {
      return this->m_data[j][i];
   }

   template <typename Int>
     const typename std::enable_if< std::is_integral<Int>::value and rank == 2, ValueType >::type &
   operator()( const Int i, const Int j ) const
   {
      return this->m_data[j][i];
   }

   /// rank 1
   template <typename Int>
      typename std::enable_if< std::is_integral<Int>::value and rank == 1, ValueType >::type &
   operator()( const Int i )
   {
      return this->m_data[i];
   }

   template <typename Int>
     const typename std::enable_if< std::is_integral<Int>::value and rank == 1, ValueType >::type &
   operator()( const Int i ) const
   {
      return this->m_data[i];
   }

   /// rank 1 ... to match normal access format.
   template <typename Int>
      typename std::enable_if< std::is_integral<Int>::value and rank == 1, ValueType >::type &
   operator[]( const Int i )
   {
      return this->m_data[i];
   }

   template <typename Int>
     const typename std::enable_if< std::is_integral<Int>::value and rank == 1, ValueType >::type &
   operator[]( const Int i ) const
   {
      return this->m_data[i];
   }

   template <int Dim>
   const typename std::enable_if< (Dim < rank), SliceType< StaticArrayType, Dim > >::type slice( const int index = 0 ) const
   {
      return SliceType< StaticArrayType, Dim >( *this, index );
   }

   template <int Dim>
         typename std::enable_if< (Dim < rank), SliceType< StaticArrayType, Dim > >::type slice( const int index = 0 ) 
   {
      return SliceType< StaticArrayType, Dim >( *this, index );
   }
};

template <class ArrayType, int _Dim>
struct SliceType
{
   typedef typename getBaseValueType<ArrayType>::type ValueType;

   enum { base_rank = ArrayType::rank };
   enum { rank = 1 };
   enum { dim = _Dim };

   enum { I0 = ( _Dim == 0 ) ? int(ArrayType::I0) : int(ArrayType::I1) };

   ArrayType &m_base_array;
   int m_index;

   SliceType (       ArrayType& array, const int index = 0 ) : m_index(index), m_base_array( array ) {}
   SliceType ( const ArrayType& array, const int index = 0 ) : m_index(index), m_base_array( const_cast<ArrayType&>(array) ) {}

   template <typename Int>
      typename std::enable_if< std::is_integral<Int>::value and base_rank == 2, ValueType >::type &
   operator[]( const Int i )
   {
      if ( dim == 0 )
         return this->m_base_array(i,this->m_index);
      else
         return this->m_base_array(this->m_index,i);
   }

   template <typename Int>
   const typename std::enable_if< std::is_integral<Int>::value and base_rank == 2, ValueType >::type & 
   operator[]( const Int i ) const
   {
      if ( dim == 0 )
         return this->m_base_array(i,this->m_index);
      else
         return this->m_base_array(this->m_index,i);
   }

   template <typename Int>
           typename std::enable_if< std::is_integral<Int>::value and base_rank == 2, ValueType >::type &
   operator()( const Int i )
   {
      if ( dim == 0 )
         return this->m_base_array(i,this->m_index);
      else
         return this->m_base_array(this->m_index,i);
   }

   template <typename Int>
     const typename std::enable_if< std::is_integral<Int>::value and base_rank == 2, ValueType >::type &
   operator()( const Int i ) const
   {
      if ( dim == 0 )
         return this->m_base_array(i,this->m_index);
      else
         return this->m_base_array(this->m_index,i);
   }
};

template <typename T>
struct getBaseValueType { typedef typename std::enable_if< std::is_fundamental<T>::value, T >::type type; };

template <typename T>
struct getBaseValueType< StaticArrayType<T> > { typedef typename StaticArrayType<T>::ValueType type; };

template <typename T, int Dim>
struct getBaseValueType< SliceType<T,Dim> > { typedef typename SliceType<T,Dim>::ValueType type; };


} // end namespace

#endif
