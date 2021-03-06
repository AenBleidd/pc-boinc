/**
	Copyright (c) 2016, All Rights Reserved.

	This software is in the public domain, furnished "as is", without technical
	support, and with no warranty, express or implied, as to its usefulness for
	any purpose.

	SseDoubleVectorTraits.hpp
	File with helper struct with definitions of various SSE2 SIMD operations.

	University of Trento,
	Department of Information Engineering and Computer Science

	Authors: (alphabetically ordered) Francesco Asnicar, Luca Masera,
			 Paolo Morettin, Nadir Sella, Thomas Tolio.

	Optimizations by Daniel Fruzynski
*/

#ifndef _SSE_DOUBLE_VECTOR_TRAITS_HPP_
#define _SSE_DOUBLE_VECTOR_TRAITS_HPP_

#include <emmintrin.h>
#include "VectorTypeTraitsBase.hpp"

template<typename T>
struct ScalarTypeTraits;

struct SseDoubleVectorTraits : public VectorTypeTraitsBase<SseDoubleVectorTraits>
{
	typedef double ElementType;
	typedef __m128d VectorType;
	
	typedef double HalfVectorType;
	typedef ScalarTypeTraits<ElementType> HalfVectorTraits;
	
	static const size_t ElementSize = sizeof(ElementType);
	static const size_t VectorSize = 16;
	static const size_t ElementCount = VectorSize / ElementSize;
	static const size_t DataAlignment = 16;
	
	// load vector from memory
	static VectorType load(const ElementType* mem_addr)
	{
		return _mm_load_pd(mem_addr);
	}
	
	// store vector in memory
	static void store(ElementType* mem_addr, const VectorType& a)
	{
		_mm_store_pd(mem_addr, a);
	}
	
	// create a vector with all elements equal to 0.0
	static VectorType setzero()
	{
		return _mm_setzero_pd();
	}
	
	// create a vector with all elements equal to A
	static constexpr VectorType set1(const ElementType A)
	{
		//return _mm_set1_pd(A);
		return VectorType{ A, A };
	}
	
	// create a vector with given elements
	static VectorType set(const ElementType A, const ElementType B)
	{
		return _mm_set_pd(A, B);
	}
	
	// returns a+b
	static VectorType add(const VectorType& a, const VectorType& b)
	{
		return _mm_add_pd(a, b);
	}
	
	// returns a-b
	static VectorType sub(const VectorType& a, const VectorType& b)
	{
		return _mm_sub_pd(a, b);
	}
	
	// returns a*b
	static VectorType mul(const VectorType& a, const VectorType& b)
	{
		return _mm_mul_pd(a, b);
	}
	
	// returns a/b
	static VectorType div(const VectorType& a, const VectorType& b)
	{
		return _mm_div_pd(a, b);
	}
	
	// returns sqrt(a) (square root)
	static VectorType sqrt(const VectorType& a)
	{
		return _mm_sqrt_pd(a);
	}
	
	// returns a*b+c
	static VectorType mul_add(const VectorType& a, const VectorType& b, const VectorType& c)
	{
		return add(mul(a, b), c);
	}
	
	// returns -(a*b)+c
	static VectorType neg_mul_add(const VectorType& a, const VectorType& b, const VectorType& c)
	{
		// -(a*b)+c == c - (a * b)
		return sub(c, mul(a, b));
	}
	
	static HalfVectorType get_lower_half(const VectorType& a)
	{
		return a[0];
	}
	
	static HalfVectorType get_upper_half(const VectorType& a)
	{
		return a[1];
	}
	
	// returns sum of vector elements
	static ElementType sum_elements(const VectorType& a)
	{
		return a[0] + a[1];
	}
	
	static ElementType* alloc(size_t element_count)
	{
		element_count = align_to_element_count(element_count);
		/*return static_cast<ElementType*>(aligned_alloc(
			DataAlignment, ElementSize * element_count));*/
		return static_cast<ElementType*>(_mm_malloc(
			ElementSize * element_count, DataAlignment));
	}
	
	static void free(ElementType* ptr)
	{
		//::free(ptr);
		_mm_free(ptr);
	}
	
	static ElementType get_at(const VectorType& v, size_t index)
	{
		return v[index];
	}
};

#endif // _SSE_DOUBLE_VECTOR_TRAITS_HPP_
