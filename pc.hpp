/**
    Copyright (c) 2013, All Right Reserved
    
    This software is in the public domain, furnished "as is", without technical
    support, and with no warranty, express or implied, as to its usefulness for
    any purpose.
    
    pc.hpp
    Library for the main functions of the PC algorithm.
    
    University of Trento,
    Department of Information Engineering and Computer Science
    
    Trento, fall 2013 - spring 2014

    Authors: (alphabetically ordered) Francesco Asnicar, Luca Masera,
             Paolo Morettin, Nadir Sella, Thomas Tolio.
*/

#include <utility>

#ifndef _PC
#define _PC

#define FILE_SCORE "score.txt"

// Type definition for the pair of integers.
typedef std::pair<int,int> intpair;

// Tests the causality and, when appropriate, cuts both the edges r->c and c->r.
void testAndRemove(double, Graph* &, const int, const int, const int, const double);

// Computes the correlation coefficient related to the edge r->c considering subsets of cardinality l of neaighbour edges.
double getCorrelationCoefficient(const int*, const int*, const int, Graph* &, const int, const int, double**);

// 
void iterativeComb(int*, const int, const int, Graph* &, const int, const int, const double, double**);

// Finds all the subsets adj(i)\{j} with cardinality equals to l (formally, |adj(i)\{j}| = l).
void findAllSubsets(Graph* &, const int, const int, const int, const double, int*, double**);

#endif //_PC
