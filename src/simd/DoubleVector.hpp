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
