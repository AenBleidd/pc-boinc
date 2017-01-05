/**
	Copyright (c) 2016, All Rights Reserved.

	This software is in the public domain, furnished "as is", without technical
	support, and with no warranty, express or implied, as to its usefulness for
	any purpose.

	DoubleVector.hpp
	This file checks which SIMD instructions are available, and defines various
	types and defines which will allows to best utilize available instruction
	sets.

	University of Trento,
	Department of Information Engineering and Computer Science

	Authors: (alphabetically ordered) Francesco Asnicar, Luca Masera,
			 Paolo Morettin, Nadir Sella, Thomas Tolio.

	Optimizations by Daniel Fruzynski
*/

#ifndef _DOUBLE_VECTOR_HPP_
#define _DOUBLE_VECTOR_HPP_

#ifdef __AVX__
#include "AvxDoubleVectorTraits.hpp"
#endif
#ifdef __SSE2__
#include "SseDoubleVectorTraits.hpp"
#endif
#include "ScalarTypeTraits.hpp"
#include "Vector.hpp"


#ifdef __AVX__
typedef Vector<AvxDoubleVectorTraits> AvxDoubleVector;
#endif
#ifdef __SSE2__
typedef Vector<SseDoubleVectorTraits> SseDoubleVector;
#endif
typedef Vector<ScalarTypeTraits<double>> ScalarDoubleVector;

#ifdef __AVX__
typedef AvxDoubleVector DoubleVector4;
#define HAS_DOUBLE_VECTOR_4 1
#endif
#ifdef __SSE2__
typedef SseDoubleVector DoubleVector2;
#define HAS_DOUBLE_VECTOR_2 1
#endif
typedef ScalarDoubleVector DoubleVector1;

#ifdef __AVX__
typedef AvxDoubleVector DoubleVector;
#elif defined (__SSE2__)
typedef SseDoubleVector DoubleVector;
#else
typedef ScalarDoubleVector DoubleVector;
#endif

#endif // _DOUBLE_VECTOR_HPP_
